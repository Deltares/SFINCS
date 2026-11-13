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
- Your employer or client owns or may claim rights in what you write (this is common for staff at universities, companies and government agencies): your organisation signs the [Entity CLA](.github/cla/entity_cla.md) first. Ask your manager or legal department, then send the signed agreement to sfincs@deltares.nl. Your GitHub account is then added to the agreement.

**How to sign.** When you open a pull request, the CLA bot posts a comment. Read the agreement it links to and reply on the pull request with exactly:

```
I have read the SFINCS Individual CLA and I hereby sign it.
```

The check turns green and the pull request can be reviewed and merged. Every commit author must have signed, so make sure the e-mail address in your commits is [linked to your GitHub account](https://docs.github.com/en/account-and-profile/setting-up-and-managing-your-personal-account-on-github/managing-email-preferences/adding-an-email-address-to-your-github-account).

Deltares employees do not need to sign.

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
