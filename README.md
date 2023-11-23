# Clustintime, a computational and visualization tool for time clustering of fMRI data.

<p align="center">
<img src=/logo.png height=250 width=250/>
</p>

Welcome! This toolbox contains all the functions required to run the whole time clustering pipeline from 
one command.


## Motivation

In certain experimental contexts, particularly in naturalistic paradigms, clinical experiments (such as epilepsy), 
and during resting states, the precise timing of neuronal activity remains elusive. 
To address this challenge, researchers often turn to clustering algorithms in the time domain. 
This approach enables the identification of spatial patterns in neural activity by grouping neurons or brain 
regions together based on the similarity of their temporal dynamics.

The concept of similarity in time domain is pivotal to this methodology. It involves assessing the degree of
resemblance or corresponding between temporal patterns of different signals or data points.
When we talk about measuring similarity in the time domain, we are investigating how closely the temporal profiles
of 2 or more signals align over a specified period.

By employing clustering algorithms sensitive to temporal dynamics, researchers can uncover 
meaningful patterns in neural activity, shedding light on how different neurons or brain regions organize their 
responses over time. 


### Key points of time domain clustering algorithms

1. Flexibility. They adapt to the inherent variability in timing, allowing the identification of patterns without
strict temporal contraints.
2. Pattern discovery that may have been overlooked in time-locked analysis.
3. Applicability across experimental designs, from exploratory research of spontaneous neural dynamics to controlled
tasks.


## How to install your toolbox



### Software and toolbox version requirements

These are the library version clustintime toolbox uses: 

``` python
certifi == 2022.9.24
charset-normalizer == 2.1.1
citeproc-py == 0.6.0
duecredit == 0.9.1
idna == 3.4
joblib == 1.2.0
lxml == 4.9.1
nibabel == 4.0.2
nilearn == 0.9.2
numpy == 1.23.4
packaging == 21.3
pandas == 1.5.1
pyparsing == 3.0.9
python-dateutil == 2.8.2
pytz == 2022.5
requests == 2.28.1
scikit-learn == 1.1.2
scipy == 1.9.3
six == 1.16.0
threadpoolctl == 3.1.0
urllib3 == 1.26.12
```

> [!TIP]
> In this repository, there is a requirements.txt file which contains the list above this line. In order to install all the libraries in a single command
line, execute:

``` bash
pip install -r requirements.txt
```

## Example of basic usage

Once you have installed the toolbox, 


## Getting involved


## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="20%"><a href="https://github.com/Cristina-Tobias"><img src="https://github.com/Cristina-Tobias.png?size=200" width="100px;" alt="Cristina Tobias"/><br /><sub><b>Cristina Tobias</b></sub></a><br /><a href="https://github.com/Cristina-Tobias/clustintime/tree/main/clustintime" title="Code">💻</a> <a title="Ideas">🤔</a> </td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/lauradefrutos"><img src="https://github.com/lauradefrutos.png?size=200" width="100px;" alt="Laura de Frutos-Sagastuy"/><br /><sub><b>Laura de Frutos</b></sub></a><br /><a href="https://github.com/ME-ICA/tedana/commits?author=dowdlelt" title="Code">💻</a> <a href="https://github.com/lauradefrutos/clustintime/blob/Documentation/README.md" title="Documentation">💬</a> </td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/lmansoo"><img src="https://github.com/lmansoo.png?size=200" width="100px;" alt="Lucia Manso-Ortega"/><br /><sub><b>Lucia Manso-Ortega</b></sub></a><br /><a href="https://github.com/lmansoo/clustintime" title="Code">💻</a> <a title="Testing">⚠️</a> </td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/marco7877"><img src="https://avatars.githubusercontent.com/u/56403434?v=4?s=100" width="100px;" alt="Marco Flores-Coronado"/><br /><sub><b>Marco Flores-Coronado</b></sub></a><br /><a href="#ideas-marco7877" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/ME-ICA/tedana/commits?author=marco7877" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/javiergcas"><img src="https://avatars1.githubusercontent.com/u/7314358?v=4?s=100" width="100px;" alt="Javier Gonzalez-Castillo"/><br /><sub><b>Javier Gonzalez-Castillo</b></sub></a><br /><a href="#ideas-javiergcas" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/ME-ICA/tedana/commits?author=javiergcas" title="Code">💻</a> <a href="#design-javiergcas" title="Design">🎨</a></td>
      <td align="center" valign="top" width="20%"><a href="https://github.com/handwerkerd"><img src="https://avatars3.githubusercontent.com/u/7406227?v=4?s=100" width="100px;" alt="Dan Handwerker"/><br /><sub><b>Dan Handwerker</b></sub></a><br /><a href="#design-handwerkerd" title="Design">🎨</a> <a href="https://github.com/ME-ICA/tedana/commits?author=handwerkerd" title="Documentation">📖</a> <a href="#example-handwerkerd" title="Examples">💡</a> <a href="https://github.com/ME-ICA/tedana/pulls?q=is%3Apr+reviewed-by%3Ahandwerkerd" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/ME-ICA/tedana/commits?author=handwerkerd" title="Code">💻</a></td>
    </tr>    
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

## License
GNU GENERAL PUBLIC LICENSE (Version 3, 29 June 2007)

Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
Everyone is permitted to copy and distribute verbatim copies
of this license document, but changing it is not allowed.
