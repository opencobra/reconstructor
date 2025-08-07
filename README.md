<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->

<a id="readme-top"></a>

<!--
*** This README was generated from the "Best‑README‑Template" with project‑specific details for the
*** openCOBRA *Constructor* application.
*** Feel free to keep iterating – just ask and we can refine any part! :)
-->

<!-- PROJECT SHIELDS -->

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![GPL‑3.0 License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- PROJECT LOGO -->

<br />
<div align="center">

  <h3 align="center">Constructor – Reconstruction Curation Tool</h3>

  <p align="center">
    A web application to build, curate, and validate biochemical reactions for genome‑scale reconstructions.
    <br />
    <a href="https://constructor.humanmetabolism.org"><strong>Launch Constructor »</strong></a>
    <br />
    <br />
    <a href="https://github.com/opencobra/reconstructor">View Source</a>
    &nbsp;·&nbsp;
    <a href="https://github.com/opencobra/reconstructor/issues/new?labels=bug&template=bug_report.md">Report Bug</a>
    &nbsp;·&nbsp;
    <a href="https://github.com/opencobra/reconstructor/issues/new?labels=enhancement&template=feature_request.md">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#key-features">Key Features</a></li>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li><a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
    <li><a href="#citation">Citation</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->

## About The Project

[![Constructor Screenshot][product-screenshot]](https://constructor.humanmetabolism.org)

Constructor (sometimes referred to simply as the _Reconstruction Curation Tool_) helps domain experts and the wider metabolic‑modeling community create **balanced, well‑annotated biochemical reactions** that can be seamlessly incorporated into genome‑scale metabolic reconstructions such as those hosted by the [Virtual Metabolic Human (VMH)](https://vmh.life/).

The traditional workflow of curating reactions across multiple spreadsheets and scripts is **time‑consuming and error‑prone**. Constructor brings everything into a single, collaborative interface:

- Manage substrates/products, compartments, charge and atom balance in real‑time.
- Fetch existing reactions from **VMH** and **Rhea**, or start from scratch.
- Draw metabolites with **ChemDoodle**, or paste identifiers (ChEBI, VMH, etc.).
- Automatically assess mass & charge balance and highlight discrepancies.
- Annotate reactions with references, external links, gene‑protein‑reaction (GPR) rules and organ localisation.
- Save reactions privately, share with collaborators, or **push directly to VMH** with one click.
- Keep track of community activity with a built‑in leaderboard and statistics dashboard.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Key Features

- 🧪 **Interactive Reaction Builder** – drag‑and‑drop UI with on‑the‑fly validation.
- 🔗 **Database Integration** – import from VMH/Rhea and export back to VMH.
- 🧬 **GPR & Gene‑Info Parsing** – associate genes, organs and sub‑cellular localisation.
- 📊 **Leaderboard & Stats** – Chart.js powered dashboard of community curation.
- 🖼️ **ChemDoodle & 3Dmol.js** – visualise 2‑D sketches and 3‑D structures inline.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With

This stack aims to be familiar to both Python/Django developers and front‑end contributors.

- [![Python][Python]][Python-url]
- [![Django][Django]][Django-url]
- ![PostgreSQL][Postgres]
- ![jQuery][JQuery.com]
- ![Semantic UI][Semantic]
- ![Chart.js][Chartjs]
- ![ChemDoodle][ChemDoodle]
- ![3Dmol.js][ThreeDmol]

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->

## Getting Started

Follow these instructions to get a local development instance running.

### Prerequisites

- **Python >= 3.10**
- **Django >= 5.0**
- **PostgreSQL** (or another Django‑compatible RDBMS)
- **Reaction Decoder Tool (RDT) >= 2.4.1** – [Download](https://github.com/asad/ReactionDecoder/releases), and put the jar file under folder curationTool
- **Matlab 2024a**

```bash
# Ubuntu example – install system packages
sudo apt update && sudo apt install python3.10 python3.10-venv build-essential postgresql postgresql-contrib
```

### Installation

```bash
# 1. Clone the repo
$ git clone https://github.com/opencobra/reconstructor.git && cd reconstructor

# 2. Create and activate a virtual environment
$ python3.10 -m venv .venv
$ source .venv/bin/activate

# 3. Install Python dependencies
$ pip install -r requirements.txt

# 4. Configure database and create config.json in root directory
# Follow the format in example.config.json, and fill in requested values

# 5. Apply migrations and run
$ cd curationTool
$ python manage.py makemigrations
$ python manage.py migrate
$ python manage.py runserver
```

Visit [http://127.0.0.1:8000/](http://127.0.0.1:8000/) to start curating reactions 🌱.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->

## Usage

Below is the most common workflow – see the [tutorials](./tutorials) directory for step‑by‑step guides.

1. **Create or fetch** a reaction in the _Reactants_ tab.
2. Inspect **atom mapping** and **chemical balance** alerts.
3. Add **references**, **external links**, **comments** or **gene rules** in their respective tabs.
4. **Save** the reaction (or _Save As_ to duplicate).
5. Navigate to _Saved Reactions_ → select entries → **Add to VMH**.

For code structure explanation, please check [docs/CODE_STRUCTURE.md](./tutorials/docs/CODE_STRUCTURE.md).

_For more examples, please refer to the [Documentation](https://opencobra.github.io)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->

## Contributing

Contributions keep Constructor improving – whether it's code, documentation, or testing. All PRs are welcome!

1. **Fork** the project
2. Create a feature branch `git checkout -b feature/awesome‑feature`
3. Commit your changes `git commit -m 'feat: add awesome feature'`
4. Push to your branch `git push origin feature/awesome‑feature`
5. Open a **pull request**

### Community leaderboard

<a href="https://github.com/opencobra/reconstructor/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=opencobra/reconstructor" alt="contributors"/>
</a>

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->

## License

Distributed under the **GNU General Public License v3.0**. See `LICENSE` for details.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->

## Contact

openCOBRA Team – [cobra@sysbio.org](mailto:cobra@sysbio.org)

DMTC Team - [digitalmetabolictwin.org](mailto:ronan.mt.fleming@universityofgalway.ie)

Project Link: [https://github.com/opencobra/reconstructor](https://github.com/opencobra/reconstructor)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ACKNOWLEDGMENTS -->

## Acknowledgments

- [Virtual Metabolic Human](https://vmh.life/)
- [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder)
- [ChemDoodle Web Components](https://web.chemdoodle.com/) and [3Dmol.js](https://3dmol.org/)
- [Chart.js](https://www.chartjs.org/)
- [Semantic UI](https://semantic-ui.com/)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Citation

Cite this work in your paper:

Alwer, S., Sathe, V., Brennan, A., and McGoldrick J. "Constructor: An open-source interface for quality-controlled metabolic reconstruction." Available at: https://github.com/opencobra/reconstructor.

<!-- MARKDOWN LINKS & IMAGES -->

<!-- Shields -->

[contributors-shield]: https://img.shields.io/github/contributors/opencobra/reconstructor.svg?style=for-the-badge
[contributors-url]: https://github.com/opencobra/reconstructor/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/opencobra/reconstructor.svg?style=for-the-badge
[forks-url]: https://github.com/opencobra/reconstructor/network/members
[stars-shield]: https://img.shields.io/github/stars/opencobra/reconstructor.svg?style=for-the-badge
[stars-url]: https://github.com/opencobra/reconstructor/stargazers
[issues-shield]: https://img.shields.io/github/issues/opencobra/reconstructor.svg?style=for-the-badge
[issues-url]: https://github.com/opencobra/reconstructor/issues
[license-shield]: https://img.shields.io/github/license/opencobra/reconstructor.svg?style=for-the-badge
[license-url]: https://github.com/opencobra/reconstructor/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/company/virtual-metabolic-human

<!-- Stack badges -->

[Python]: https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white
[Python-url]: https://www.python.org/
[Django]: https://img.shields.io/badge/Django-092E20?style=for-the-badge&logo=django&logoColor=white
[Django-url]: https://www.djangoproject.com/
[Postgres]: https://img.shields.io/badge/PostgreSQL-4169E1?style=for-the-badge&logo=postgresql&logoColor=white
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[Semantic]: https://img.shields.io/badge/Semantic%20UI-35BDB2?style=for-the-badge&logo=semanticui&logoColor=white
[Chartjs]: https://img.shields.io/badge/Chart.js-F5788D?style=for-the-badge&logo=chartdotjs&logoColor=white
[ChemDoodle]: https://img.shields.io/badge/ChemDoodle-2780D4?style=for-the-badge
[ThreeDmol]: https://img.shields.io/badge/3Dmol.js-EB4E00?style=for-the-badge

<!-- Images -->

[product-screenshot]: Reconstruction-Interface-Architecture.png
