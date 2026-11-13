# Contributing to SFINCS

Thank you for your interest in SFINCS. Contributions of all kinds are welcome: bug reports, documentation, test cases and code.

SFINCS is owned and maintained by [Deltares](https://www.deltares.nl). Please read this page before opening a pull request.

## Reporting issues and asking questions

- Bugs and feature requests: the [issue tracker](https://github.com/Deltares/SFINCS/issues).
- Questions about using SFINCS: the [discussion board](https://github.com/Deltares/SFINCS/discussions).
- Questions about building models with HydroMT-SFINCS: the [HydroMT-SFINCS discussion board](https://github.com/Deltares/hydromt_sfincs/discussions).
- Anything else: sfincs@deltares.nl.

Reporting an issue or suggesting an idea does not require a contributor agreement.

## Contributor License Agreement (CLA)

Before we can merge your first pull request, you need to accept the SFINCS Contributor License Agreement. You only do this once; it then covers all your contributions.

**Why?** You keep the copyright in your work. The CLA gives Deltares a permanent licence to use and distribute your contribution, including under licences other than the current GPL-3.0. Deltares already distributes SFINCS under more than one licence (the GPL-3.0 source code and the freeware executables), and needs this right to keep doing so.

**Which agreement?**

- You contribute in your own time, or your employer has no rights in your work: the [Individual CLA](.github/cla/individual_cla.md).
- Your employer or client owns or may claim rights in what you write (this is common for staff at universities, companies and government agencies): your organisation signs the [Entity CLA](.github/cla/entity_cla.md) first. Ask your manager or legal department, then send the signed agreement to sfincs@deltares.nl. Unless Deltares tells you otherwise, you then do not sign the individual CLA; your GitHub account is added to the exemption list instead (see below).

**How to sign.** When you open a pull request, the CLA bot posts a comment. Read the agreement it links to and reply on the pull request with exactly:

```
I have read the SFINCS Individual CLA and I hereby sign it.
```

The check turns green and the pull request can be reviewed and merged. Every commit author must have signed, so make sure the e-mail address in your commits is [linked to your GitHub account](https://docs.github.com/en/account-and-profile/setting-up-and-managing-your-personal-account-on-github/managing-email-preferences/adding-an-email-address-to-your-github-account).

### Exemptions

The CLA bot does not ask the following GitHub accounts to sign. They are listed in the `allowlist` of [.github/workflows/cla.yml](.github/workflows/cla.yml):

- **Deltares employees.** Work you do for Deltares already belongs to Deltares, so you do not sign the CLA. If the bot asks you to sign, do not sign; ask a maintainer to add your GitHub user name to the list.
- **Designated contributors of an organisation that signed the entity CLA.** When Deltares receives a signed entity CLA, a maintainer adds the GitHub accounts in its Schedule A to the list. When the organisation reports that someone has left, the maintainer removes the account.
- **Automation accounts** such as Dependabot and Copilot.

Changes to the list go through a normal pull request reviewed by a Deltares maintainer. People who leave Deltares are removed from the list; their later contributions need a CLA like anyone else's.

## Code from elsewhere

Only submit work you wrote yourself. If a contribution includes code from another project (a library, a published algorithm implementation, code generated from someone else's source), say so in the pull request, with its origin and licence. Do not add your own copyright lines to source file headers; the Deltares copyright notice stays in place.

## Pull requests

- Branch from `main` and keep each pull request focused on one change.
- Fill in the checklist in the pull request template.
- Fortran changes follow the [coding style](source/CODING_STYLE.md).
- Changes to model functionality are tested with the SFINCS testbed; new functionality comes with a conceptual and a real-world test case.
- Update the [documentation](docs) when you change input, output or behaviour.

## Licence

SFINCS source code is distributed under the [GNU GPL-3.0](LICENSE). Contributions accepted under the CLA may also be distributed by Deltares under other licences, as described above.
