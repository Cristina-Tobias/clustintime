# Clustintime, a computational and visualization tool for time clustering of fMRI data.

Welcome! This toolbox contains all the functions required to run the whole time clustering pipeline from 
one command.

<img src=/logo.png height=100 width=100/>


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


> [!TIP]
> In this repository, there is a requirements.txt file which contains the list above this line. In order to install all the libraries in a single command
line, execute:

``` ruby
pip install -r requirements.txt
```

## Example of basic usage

Once you have installed the toolbox, 


## Include badges of build status, code coverage and any other relevant metrics.


## License

