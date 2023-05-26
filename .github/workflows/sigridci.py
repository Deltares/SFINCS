#!/usr/bin/env python3

# Copyright Software Improvement Group
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import base64
import datetime
import dataclasses
import html
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import typing
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass
from enum import Enum
from xml.dom import minidom


LOG_HISTORY = []


def log(message):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{timestamp}  {message}", flush=True)
    LOG_HISTORY.append(message)


@dataclass
class UploadOptions:
    sourceDir: str = None
    excludePatterns: typing.List[str] = dataclasses.field(default_factory=lambda: [])
    includeHistory: bool = False
    subsystem: str = ""
    showContents: bool = False
    publishOnly: bool = False
    
    def readScopeFile(self):
        return self.locateFile(["sigrid.yaml", "sigrid.yml"])
        
    def readMetadataFile(self):
        return self.locateFile(["sigrid-metadata.yaml", "sigrid-metadata.yml"])
    
    def locateFile(self, possibleFileNames):
        for file in possibleFileNames:
            if os.path.exists(f"{self.sourceDir}/{file}"):
                with open(f"{self.sourceDir}/{file}", "r") as f:
                    return f.read()
        return None


class TargetQuality:
    def __init__(self, scope, targetRating):
        self.ratings = {"MAINTAINABILITY" : targetRating}
        
        # We can't use pyyaml because PIP is not available in the some of the
        # very diverse set of customer environments where Sigrid CI is used.
        targetPattern = re.compile("(" + "|".join(Report.METRICS) + "):\s*([\d\.]+)", re.IGNORECASE)
        
        for line in scope.split("\n"):
            match = targetPattern.match(line.strip())
            if match:
                log(f"Loading {match.group(1).upper()} target from scope configuration file")
                self.ratings[match.group(1).upper()] = float(match.group(2))

    def meetsTargetQualityForMetric(self, feedback, metric):
        value = feedback["newCodeRatings"].get(metric, None)
        targetRating = self.ratings.get(metric, None)
        return value == None or targetRating == None or value >= targetRating
        
    def meetsQualityTargets(self, feedback):
        return all(self.meetsTargetQualityForMetric(feedback, metric) for metric in self.ratings)
        
        
class RunMode(Enum):
    FEEDBACK_ONLY = 1
    FEEDBACK_AND_PUBLISH = 2
    PUBLISH_ONLY = 3


class SigridApiClient:
    API_VERSION = "v1"
    POLL_INTERVAL = 30
    POLL_ATTEMPTS = 120

    def __init__(self, args):
        self.baseURL = args.sigridurl
        self.urlPartnerName = urllib.parse.quote_plus(args.partner.lower())
        self.urlCustomerName = urllib.parse.quote_plus(args.customer.lower())
        self.urlSystemName = urllib.parse.quote_plus(args.system.lower())
        self.subsystem = args.subsystem
        
        if args.publishonly:
            self.runMode = RunMode.PUBLISH_ONLY
        elif args.publish:
            self.runMode = RunMode.FEEDBACK_AND_PUBLISH
        else:
            self.runMode = RunMode.FEEDBACK_ONLY

    def callSigridAPI(self, path, body=None, contentType=None):
        url = f"{self.baseURL}/rest/{path}"
        request = urllib.request.Request(url, body)
        request.add_header("Accept", "application/json")
        request.add_header("Authorization", f"Bearer {os.environ['SIGRID_CI_TOKEN']}".encode("utf8"))
        if contentType != None:
            request.add_header("Content-Type", contentType)

        response = urllib.request.urlopen(request)
        if response.status == 204:
            return {}
        responseBody = response.read().decode("utf8")
        if len(responseBody) == 0:
            log("Received empty response")
            return {}
        return json.loads(responseBody)
        
    def retry(self, operation, *, attempts=5, allow404=False, allowEmpty=True):
        for attempt in range(attempts):
            try:
                response = operation()
                if allowEmpty or response != {}:
                    return response
            except urllib.error.HTTPError as e:
                if e.code in [401, 403]:
                    log("You are not authorized to access Sigrid for this system")
                    sys.exit(1)
                elif allow404 and e.code == 404:
                    return False
            
            # These statements are intentionally outside of the except-block,
            # since we want to retry for empty response on some end points.
            log("Retrying")
            time.sleep(self.POLL_INTERVAL)
        
        log(f"Sigrid is currently unavailable, failed after {attempts} attempts")
        sys.exit(1)

    def submitUpload(self, options, systemExists):
        with tempfile.TemporaryDirectory() as tempDir:
            log("Creating upload")
            uploadPacker = SystemUploadPacker(options)
            upload = os.path.join(tempDir, "upload.zip")
            uploadPacker.prepareUpload(options.sourceDir, upload)

            log("Preparing upload")
            uploadLocation = self.obtainUploadLocation(systemExists)
            uploadUrl = uploadLocation["uploadUrl"]
            analysisId = uploadLocation["ciRunId"]
            log(f"Sigrid CI analysis ID: {analysisId}")
            log("Submitting upload" if self.runMode == RunMode.FEEDBACK_ONLY else "Publishing upload")
            self.uploadBinaryFile(uploadUrl, upload)

            return analysisId

    def obtainUploadLocation(self, systemExists):
        path = f"/inboundresults/{self.urlPartnerName}/{self.urlCustomerName}/{self.urlSystemName}/ci/uploads/{self.API_VERSION}"
        if not systemExists:
            path += "/onboarding"
        elif self.runMode == RunMode.PUBLISH_ONLY:
            path += "/publishonly"
        elif self.runMode == RunMode.FEEDBACK_AND_PUBLISH:
            path += "/publish"
        
        if self.subsystem:
            path += "?subsystem=" + urllib.parse.quote_plus(self.subsystem)
    
        return self.retry(lambda: self.callSigridAPI(path))
        
    def validateScopeFile(self, scopeFile):
        path = f"/inboundresults/{self.urlPartnerName}/{self.urlCustomerName}/{self.urlSystemName}/ci/validate/{self.API_VERSION}"
        return self.retry(lambda: self.callSigridAPI(path, scopeFile.encode("utf8"), "application/yaml"))
        
    def validateMetadata(self, metadataFile):
        path = f"/analysis-results/sigridci/{self.urlCustomerName}/validate"
        return self.retry(lambda: self.callSigridAPI(path, metadataFile.encode("utf8"), "application/yaml"))

    def uploadBinaryFile(self, url, upload):
        self.retry(lambda: self.attemptUpload(url, upload))
        log(f"Upload successful")

    def attemptUpload(self, url, upload):
        with open(upload, "rb") as uploadRef:
            uploadRequest = urllib.request.Request(url, data=uploadRef)
            uploadRequest.method = "PUT"
            uploadRequest.add_header("Content-Type", "application/zip")
            uploadRequest.add_header("Content-Length", "%d" % os.path.getsize(upload))
            uploadRequest.add_header("x-amz-server-side-encryption", "AES256")
            urllib.request.urlopen(uploadRequest)

    def checkSystemExists(self):
        path = f"/analysis-results/sigridci/{self.urlCustomerName}/{self.urlSystemName}/{self.API_VERSION}/ci"
        return self.retry(lambda: self.callSigridAPI(path), allow404=True) != False

    def fetchAnalysisResults(self, analysisId):
        log("Waiting for analysis results")
        path = f"/analysis-results/sigridci/{self.urlCustomerName}/{self.urlSystemName}/{self.API_VERSION}/ci/results/{analysisId}"
        return self.retry(lambda: self.callSigridAPI(path), attempts=self.POLL_ATTEMPTS, allowEmpty=False)
        
    def fetchMetadata(self):
        path = f"/analysis-results/api/{self.API_VERSION}/system-metadata/{self.urlCustomerName}/{self.urlSystemName}"
        return self.retry(lambda: self.callSigridAPI(path))
        
    def fetchObjectives(self):
        path = f"/analysis-results/api/{self.API_VERSION}/objectives/{self.urlCustomerName}/{self.urlSystemName}/config"
        return self.retry(lambda: self.callSigridAPI(path))
        
    def getLandingPage(self, analysisId, target):
        targetRating = "%.1f" % target.ratings["MAINTAINABILITY"]
        return f"{self.baseURL}/{self.urlCustomerName}/{self.urlSystemName}/-/sigrid-ci/{analysisId}?targetRating={targetRating}"


class SystemUploadPacker:
    MAX_UPLOAD_SIZE_MB = 500

    DEFAULT_EXCLUDES = [
        "$tf/",
        "coverage/",
        "build/",
        "dist/",
        "node_modules/",
        "sigridci/",
        "sigrid-ci-output/",
        "target/",
        ".git/",
        ".gitattributes",
        ".gitignore",
        ".idea/",
        ".jpg",
        ".png"
    ]

    def __init__(self, options):
        self.options = options

    def prepareUpload(self, sourceDir, outputFile):
        zipFile = zipfile.ZipFile(outputFile, "w", zipfile.ZIP_DEFLATED)
        hasContents = False
        
        if self.options.includeHistory and os.path.exists(f"{sourceDir}/.git"):
            self.includeRepositoryHistory(sourceDir)

        for root, dirs, files in os.walk(sourceDir):
            for file in sorted(files):
                filePath = os.path.join(root, file)
                if file != outputFile and not self.isExcluded(filePath):
                    relativePath = os.path.relpath(os.path.join(root, file), sourceDir)
                    hasContents = True
                    if self.options.showContents:
                        log(f"Adding file to upload: {relativePath}")
                    zipFile.write(filePath, relativePath)

        zipFile.close()

        self.checkUploadContents(outputFile, hasContents)

    def checkUploadContents(self, outputFile, hasContents):
        uploadSizeBytes = os.path.getsize(outputFile)
        uploadSizeMB = max(round(uploadSizeBytes / 1024 / 1024), 1)
        log(f"Upload size is {uploadSizeMB} MB")

        if uploadSizeMB > self.MAX_UPLOAD_SIZE_MB:
            raise Exception(f"Upload exceeds maximum size of {self.MAX_UPLOAD_SIZE_MB} MB")
        elif not hasContents:
            print(f"No code found to upload, please check the directory used for --source")
            sys.exit(1)
        elif uploadSizeBytes < 50000:
            log("Warning: Upload is very small, source directory might not contain all source code")

    def isExcluded(self, filePath):
        excludePatterns = self.DEFAULT_EXCLUDES + (self.options.excludePatterns or [])
        normalizedPath = filePath.replace("\\", "/")
        for exclude in excludePatterns:
            if exclude != "" and exclude.strip() in normalizedPath:
                return True
        return False
        
    def includeRepositoryHistory(self, sourceDir):
        gitCommand = ["git", "-C", sourceDir, "--no-pager", "log", "--date=iso", "--format='@@@;%H;%an;%ae;%ad;%s'", \
                      "--numstat", "--no-merges"]
        try:
            output = subprocess.run(gitCommand, stdout=subprocess.PIPE)
            if output.returncode == 0:
                with open(f"{sourceDir}/git.log", "w") as f:
                    f.write(output.stdout.decode("utf8"))
            else:
                log("Exporting repository history failed")
        except Exception as e:
            log("Error while trying to include repository history: " + str(e))
    

class Report:
    METRICS = ["VOLUME", "DUPLICATION", "UNIT_SIZE", "UNIT_COMPLEXITY", "UNIT_INTERFACING", "MODULE_COUPLING",
               "COMPONENT_BALANCE_PROP", "COMPONENT_INDEPENDENCE", "COMPONENT_ENTANGLEMENT", "MAINTAINABILITY"]

    REFACTORING_CANDIDATE_METRICS = ["DUPLICATION", "UNIT_SIZE", "UNIT_COMPLEXITY", "UNIT_INTERFACING",
                                     "MODULE_COUPLING"]

    def generate(self, analysisId, feedback, args, target):
        pass

    def formatMetricName(self, metric):
        return metric.replace("_PROP", "").title().replace("_", " ")

    def formatRating(self, ratings, metric, naText="N/A"):
        if ratings.get(metric, None) == None:
            return naText
        return "%.1f" % ratings[metric]

    def formatBaseline(self, feedback):
        if not feedback.get("baseline", None):
            return "N/A"
        snapshotDate = datetime.datetime.strptime(feedback["baseline"], "%Y%m%d")
        return snapshotDate.strftime("%Y-%m-%d")

    def getSigridUrl(self, args):
        customer = urllib.parse.quote_plus(args.customer.lower())
        system = urllib.parse.quote_plus(args.system.lower())
        return f"https://sigrid-says.com/{customer}/{system}"

    def getRefactoringCandidates(self, feedback, metric):
        refactoringCandidates = feedback.get("refactoringCandidates", [])
        return [rc for rc in refactoringCandidates if rc["metric"] == metric or metric == "MAINTAINABILITY"]


class TextReport(Report):
    ANSI_BOLD = "\033[1m"
    ANSI_GREEN = "\033[92m"
    ANSI_YELLOW = "\033[33m"
    ANSI_RED = "\033[91m"
    ANSI_BLUE = "\033[96m"
    LINE_WIDTH = 89

    def __init__(self, output=sys.stdout):
        self.output = output

    def generate(self, analysisId, feedback, args, target):
        self.printHeader("Refactoring candidates")
        for metric in self.REFACTORING_CANDIDATE_METRICS:
            self.printMetric(feedback, metric)

        self.printHeader("Maintainability ratings")
        self.printTableRow(["System property", f"Baseline on {self.formatBaseline(feedback)}", \
            "New/changed code", "Target", "Overall" if args.publish else ""] )

        for metric in self.METRICS:
            if metric == "MAINTAINABILITY":
                self.printSeparator()

            row = [
                self.formatMetricName(metric),
                "(" + self.formatRating(feedback["baselineRatings"], metric) + ")",
                self.formatRating(feedback["newCodeRatings"], metric),
                str(target.ratings.get(metric, "")),
                self.formatRating(feedback["baselineRatings"], metric) if args.publish else ""
            ]

            self.printTableRow(row, self.getRatingColor(feedback, target, metric))

    def printTableRow(self, row, color=None):
        formattedRow = "%-27s%-25s%-20s%-10s%-7s" % tuple(row)
        if color:
            self.printColor(formattedRow, color)
        else:
            print(formattedRow, file=self.output)

    def printHeader(self, header):
        print("", file=self.output)
        self.printSeparator()
        print(header, file=self.output)
        self.printSeparator()

    def printSeparator(self):
        print("-" * self.LINE_WIDTH, file=self.output)

    def printMetric(self, feedback, metric):
        print("", file=self.output)
        print(self.formatMetricName(metric), file=self.output)

        refactoringCandidates = self.getRefactoringCandidates(feedback, metric)
        if len(refactoringCandidates) == 0:
            print("    None", file=self.output)
        else:
            for rc in refactoringCandidates:
                print(self.formatRefactoringCandidate(rc), file=self.output)

    def getRatingColor(self, feedback, target, metric):
        if feedback["newCodeRatings"].get(metric, None) == None or not metric in target.ratings:
            return self.ANSI_BLUE
        elif target.meetsTargetQualityForMetric(feedback, metric):
            return self.ANSI_GREEN
        else:
            return self.ANSI_RED

    def formatRefactoringCandidate(self, rc):
        category = ("(" + rc["category"] + ")").ljust(14)
        subject = rc["subject"].replace("\n", "\n" + (" " * 21)).replace("::", "\n" + (" " * 21))
        return f"    - {category} {subject}"

    def printColor(self, message, ansiPrefix):
        print(ansiPrefix + message + "\033[0m", file=self.output)


class StaticHtmlReport(Report):
    HTML_STAR_FULL = "&#9733;"
    HTML_STAR_EMPTY = "&#9734;"

    def generate(self, analysisId, feedback, args, target):
        with open(os.path.dirname(__file__) + "/sigridci-feedback-template.html", encoding="utf-8", mode="r") as templateRef:
            template = templateRef.read()
            template = self.renderHtmlFeedback(template, feedback, args, target)

        reportFile = os.path.abspath("sigrid-ci-output/index.html")
        writer = open(reportFile, encoding="utf-8", mode="w")
        writer.write(template)
        writer.close()

        print("")
        print("You can find the full results here:")
        print("    " + reportFile)
        print("")
        print("You can find more information about these results in Sigrid:")
        print("    " + self.getSigridUrl(args))
        print("")

    def renderHtmlFeedback(self, template, feedback, args, target):
        placeholders = {
            "CUSTOMER" : html.escape(args.customer),
            "SYSTEM" : html.escape(args.system),
            "TARGET" : "%.1f" % target.ratings["MAINTAINABILITY"],
            "LINES_OF_CODE_TOUCHED" : "%d" % feedback.get("newCodeLinesOfCode", 0),
            "BASELINE_DATE" : self.formatBaseline(feedback),
            "SIGRID_LINK" : self.getSigridUrl(args),
            "MAINTAINABILITY_PASSED" : ("passed" if target.meetsQualityTargets(feedback) else "failed")
        }

        for metric in self.METRICS:
            placeholders[f"{metric}_OVERALL"] = self.formatRating(feedback["baselineRatings"], metric)
            placeholders[f"{metric}_NEW"] = self.formatRating(feedback["newCodeRatings"], metric)
            placeholders[f"{metric}_TARGET"] = self.formatRating(target.ratings, metric, "")
            placeholders[f"{metric}_STARS_OVERALL"] = self.formatHtmlStars(feedback["baselineRatings"], metric)
            placeholders[f"{metric}_STARS_NEW"] = self.formatHtmlStars(feedback["newCodeRatings"], metric)
            placeholders[f"{metric}_PASSED"] = self.formatPassed(feedback, target, metric)
            placeholders[f"{metric}_REFACTORING_CANDIDATES"] = self.formatRefactoringCandidates(feedback, metric)

        return self.fillPlaceholders(template, placeholders)

    def fillPlaceholders(self, template, placeholders):
        for placeholder, value in placeholders.items():
            template = template.replace(f"@@@{placeholder}", value)
        return template

    def formatPassed(self, feedback, target, metric):
        if target.ratings.get(metric, None) == None:
            return ""
        return "passed" if target.meetsTargetQualityForMetric(feedback, metric) else "failed"

    def formatRefactoringCandidates(self, feedback, metric):
        refactoringCandidates = self.getRefactoringCandidates(feedback, metric)
        if len(refactoringCandidates) == 0:
            return "None"
        return "\n".join([self.formatRefactoringCandidate(rc) for rc in refactoringCandidates])

    def formatRefactoringCandidate(self, rc):
        subjectName = html.escape(rc["subject"]).replace("\n", "<br />").replace("::", "<br />")
        category = html.escape(rc["category"])
        return f"<span><em>({category})</em><div>{subjectName}</div></span>"

    def formatHtmlStars(self, ratings, metric):
        if ratings.get(metric, None) == None:
            return "N/A"
        stars = min(int(ratings[metric] + 0.5), 5)
        fullStars = stars * self.HTML_STAR_FULL
        emptyStars = (5 - stars) * self.HTML_STAR_EMPTY
        rating = self.formatRating(ratings, metric)
        return f"<strong class=\"stars{stars}\">{fullStars}{emptyStars}</strong> &nbsp; " + rating


class JUnitFormatReport(Report):
    def generate(self, analysisId, feedback, args, target):
        with open("sigrid-ci-output/sigridci-junit-format-report.xml", "w") as fileRef:
            fileRef.write(self.generateXML(feedback, target))

    def generateXML(self, feedback, target):
        dom = minidom.Document()
        testSuite = dom.createElement("testsuite")
        testSuite.setAttribute("name", "Sigrid CI")
        dom.appendChild(testSuite)

        testCase = dom.createElement("testcase")
        testCase.setAttribute("classname", "Sigrid CI")
        testCase.setAttribute("name", "Maintainability")
        testSuite.appendChild(testCase)

        failures = self.getFailures(feedback, target)
        if len(failures) > 0:
            failure = dom.createElement("failure")
            failure.appendChild(dom.createTextNode("Refactoring candidates:\n\n" + "\n".join(failures)))
            testCase.appendChild(failure)

        return dom.toprettyxml(indent="    ")

    def getFailures(self, feedback, target):
        if target.meetsQualityTargets(feedback):
            return []

        formatFailure = lambda rc: f"- {rc['subject']}\n  ({self.formatMetricName(rc['metric'])}, {rc['category']})"

        candidates = self.getRefactoringCandidates(feedback, "MAINTAINABILITY")\
            if not target.meetsTargetQualityForMetric(feedback, "MAINTAINABILITY")\
            else self.getFailedRcs(feedback, target)
        return [formatFailure(rc) for rc in candidates]

    def getFailedRcs(self, feedback, target):
        lists = [self.getRefactoringCandidates(feedback, m) for m in target.ratings
                 if not target.meetsTargetQualityForMetric(feedback, m)]
        return [item for sublist in lists for item in sublist]
                
                
class ConclusionReport(Report):
    def __init__(self, apiClient, output=sys.stdout):
        self.apiClient = apiClient
        self.output = output

    def generate(self, analysisId, feedback, args, target):
        self.printConclusionMessage(feedback, target)
        self.printLandingPage(analysisId, feedback, target)
        # If you publish(only) we never break the build
        # We can break the build when running on a branch or pull request.
        if not target.meetsQualityTargets(feedback) and not args.publish:
            sys.exit(1)
            
    def printConclusionMessage(self, feedback, target):
        asciiArt = TextReport(self.output)

        if feedback["newCodeRatings"].get("MAINTAINABILITY", None) == None:
            asciiArt.printColor("\n** SIGRID CI RUN COMPLETE: NO FILES CONSIDERED FOR MAINTAINABILITY WERE CHANGED **\n", \
                asciiArt.ANSI_BOLD + asciiArt.ANSI_BLUE)
        elif target.meetsQualityTargets(feedback):
            asciiArt.printColor("\n** SIGRID CI RUN COMPLETE: YOU WROTE MAINTAINABLE CODE AND REACHED THE TARGET **\n", \
                asciiArt.ANSI_BOLD + asciiArt.ANSI_GREEN)
        else:
            asciiArt.printColor("\n** SIGRID CI RUN COMPLETE: THE CODE YOU WROTE DID NOT MEET THE TARGET FOR MAINTAINABLE CODE **\n", \
                asciiArt.ANSI_BOLD + asciiArt.ANSI_YELLOW)
                
    def printLandingPage(self, analysisId, feedback, target):
        landingPage = self.apiClient.getLandingPage(analysisId, target)
        
        print("", file=self.output)
        print("-" * (len(landingPage) + 4), file=self.output)
        print("View your analysis results in Sigrid:", file=self.output)
        print(f"    {landingPage}", file=self.output)
        print("-" * (len(landingPage) + 4), file=self.output)
        print("", file=self.output)


class SigridCiRunner:
    SYSTEM_NAME_PATTERN = re.compile("^[a-z0-9]+(-[a-z0-9]+)*$", re.IGNORECASE)
    SYSTEM_NAME_LENGTH = range(2, 65)
    METADATA_FIELDS = [
        "displayName",
        "divisionName",
        "teamNames",
        "supplierNames",
        "lifecyclePhase",
        "inProductionSince",
        "businessCriticality",
        "targetIndustry",
        "deploymentType",
        "applicationType",
        "externalID",
        "isDevelopmentOnly",
        "remark"
    ]
    
    def loadSigridTarget(self, apiClient):
        objectives = apiClient.fetchObjectives()
        targetRating = objectives.get("NEW_CODE_QUALITY", objectives.get("MAINTAINABILITY", 3.5))
        log("Using Sigrid for target rating (%.1f stars)" % targetRating)
        return targetRating

    def run(self, apiClient, options, target, reports):
        if os.path.exists(f"{options.sourceDir}/sigrid.yml"):
            log("Found sigrid.yml in repository. Did you mean sigrid.yaml?")
            sys.exit(1)
    
        systemExists = apiClient.checkSystemExists()
        log("Found system in Sigrid" if systemExists else "System is not yet on-boarded to Sigrid")
        
        self.prepareMetadata(options)
        self.validateConfigurationFiles(apiClient, options)
        analysisId = apiClient.submitUpload(options, systemExists)

        if not systemExists:
            log(f"System '{apiClient.urlSystemName}' is on-boarded to Sigrid, and will appear in sigrid-says.com shortly")
        elif options.publishOnly:
            log("Your project's source code has been published to Sigrid")
            self.displayMetadata(apiClient, options)
        else:
            feedback = apiClient.fetchAnalysisResults(analysisId)
            self.displayMetadata(apiClient, options)

            if not os.path.exists("sigrid-ci-output"):
                os.mkdir("sigrid-ci-output")

            for report in reports:
                report.generate(analysisId, feedback, args, target)
                
    def validateConfigurationFiles(self, apiClient, options):
        scope = options.readScopeFile()
        if scope:
            self.validateConfiguration(lambda: apiClient.validateScopeFile(scope), "scope configuration file")
        
        metadataFile = options.readMetadataFile()
        if metadataFile:
            self.validateConfiguration(lambda: apiClient.validateMetadata(metadataFile), "Sigrid metadata file")
    
    def validateConfiguration(self, validationCall, configurationName):
        log(f"Validating {configurationName}")
        validationResult = validationCall()
        if validationResult["valid"]:
            log("Validation passed")
        else:
            log("-" * 80)
            log(f"Invalid {configurationName}:")
            for note in validationResult["notes"]:
                log(f"    - {note}")
            log("-" * 80)
            sys.exit(1)
            
    def displayMetadata(self, apiClient, options):
        if options.readMetadataFile() == None:
            print("")
            print("Sigrid metadata for this system:")
            for key, value in apiClient.fetchMetadata().items():
                if value:
                    print(f"    {key}:".ljust(20) + str(value))
                    
    def prepareMetadata(self, options):
        getMetadataValue = lambda field: os.environ.get(field.lower(), "")
        metadata = {field: getMetadataValue(field) for field in self.METADATA_FIELDS if getMetadataValue(field)}
        
        if len(metadata) > 0:
            if options.readMetadataFile() != None:
                raise Exception("Cannot add metadata using environment variables if metadata YAML file is already used")
            
            with open(f"{options.sourceDir}/sigrid-metadata.yaml", "w") as writer:
                writer.write("metadata:\n")
                for name, value in metadata.items():
                    formattedValue = f"[\"{value}\"]" if name in ["teamNames", "supplierNames"] else f"\"{value}\""
                    writer.write(f"  {name}: {formattedValue}\n")
                
    def isValidSystemName(self, customerName, systemName):
        return self.SYSTEM_NAME_PATTERN.match(systemName) and \
            len(systemName) >= self.SYSTEM_NAME_LENGTH.start and \
            (len(systemName) + len(customerName) + 1) in self.SYSTEM_NAME_LENGTH


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--partner", type=str, default="sig")
    parser.add_argument("--customer", type=str)
    parser.add_argument("--system", type=str)
    parser.add_argument("--source", type=str)
    parser.add_argument("--targetquality", type=str, default="sigrid")
    parser.add_argument("--publish", action="store_true")
    parser.add_argument("--publishonly", action="store_true")
    parser.add_argument("--exclude", type=str, default="")
    parser.add_argument("--subsystem", type=str, default="")
    parser.add_argument("--showupload", action="store_true")
    parser.add_argument("--include-history", action="store_true")
    parser.add_argument("--sigridurl", type=str, default="https://sigrid-says.com")
    # Dummy argument used when passing false to boolean arguments.
    # BooleanOptionalAction would solve this, but requires Python 3.9+.
    parser.add_argument("--dummy", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.customer == None or args.system == None or args.source == None:
        parser.print_help()
        sys.exit(1)

    if sys.version_info.major == 2 or sys.version_info.minor < 7:
        print("Sigrid CI requires Python 3.7 or higher")
        sys.exit(1)

    if not "SIGRID_CI_TOKEN" in os.environ:
        print("Missing required environment variable SIGRID_CI_TOKEN")
        sys.exit(1)

    if not os.path.exists(args.source):
        print("Source code directory not found: " + args.source)
        sys.exit(1)

    log("Starting Sigrid CI")
    
    options = UploadOptions(args.source, args.exclude.split(","), args.include_history, args.subsystem, args.showupload, args.publishonly)
    apiClient = SigridApiClient(args)
    reports = [TextReport(), StaticHtmlReport(), JUnitFormatReport(), ConclusionReport(apiClient)]

    runner = SigridCiRunner()
    
    if not runner.isValidSystemName(args.customer, args.system):
        maxNameLength = runner.SYSTEM_NAME_LENGTH.stop - (len(args.customer) + 1)
        print(f"Invalid system name, system name should match '{runner.SYSTEM_NAME_PATTERN.pattern}' "
              f"and be {runner.SYSTEM_NAME_LENGTH.start} to {maxNameLength} characters long (inclusive).")
        sys.exit(1)
        
    targetRating = runner.loadSigridTarget(apiClient) if args.targetquality == "sigrid" else float(args.targetquality)
    target = TargetQuality(options.readScopeFile() or "", targetRating)
    runner.run(apiClient, options, target, reports)
