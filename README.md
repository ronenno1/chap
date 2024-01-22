![](https://in.bgu.ac.il/en/Labs/CNL/chap/images/logo.png)

# CHAP: Open Source Software for Processing and Analyzing Pupillometry Data


Pupil dilation is an effective indicator of cognitive and affective processes. While there are several eye-tracker systems in the market that provide effective solutions for pupil dilation measurement, there is a lack of tools for processing and analyzing the data provided by these systems. For this reason, we developed CHAP - open-source software written in MATLAB. This software provides a user-friendly interface (graphical user interface) for processing and analyzing pupillometry data. Our software creates uniform conventions for the pre-processing and analysis of pupillometry data, and provides a quick and easy-to-implement tool for researchers interested in pupillometry.



## Installation


The installation of CHAP includes 3 main steps:


1. Downloading CHAP from [here](https://github.com/ronenno1/chap/archive/master.zip
).
2. Extracting the ZIP file.
3. Running CHAP_installation.m file (that adds CHAP to the MATLAB path list).

Please note: 
* When working with EyeLink data files:
  * Users who use Windows 32-bit version should go to CHAP's folder and copy manually the edf2mat/edfapi.dll file to C:\Windows\System32
  * Users who use Windows 64-bit version should go to CHAP 's folder and copy manually the edf2mat/edfapi64.dll file to C:\Windows\System32
  * Users who use Mac should go to CHAP 's folder and unzip manually the 
edf2mat/edfapi.framework.zip file and copy edfapi.framework to /Library/Frameworks.

  More information about converting EyeLink data files into MATLAB can be found here: [https://github.com/uzh/edf-converter](https://github.com/uzh/edf-converter).
* When working with Eye Tribe data files (txt files given by the default Eye Tribe interface), python should be installed.


For further information or if you have any questions please do not hesitate to contact [us](mailto:Ronen.Hershman@uibk.ac.at).


## Citations

* Hershman, R., Henik, A., & Cohen, N. (2018). A novel blink detection method based on pupillometry noise. Behavior Research Methods, 50(1), 107-114. [https://doi.org/10.3758/s13428-017-1008-1](https://doi.org/10.3758/s13428-017-1008-1)
* Hershman, R., Henik, A., & Cohen, N. (2019). CHAP: Open-source software for processing and analyzing pupillometry data. Behavior Research Methods, 51(3), 1059-1074. [https://doi.org/10.3758/s13428-018-01190-1](https://doi.org/10.3758/s13428-018-01190-1).
* Shechter, A., Hershman, R., & Share, D. (2022). A Pupillometric Study of Developmental and Individual Differences in Cognitive Effort in Visual Word Recognition. Scientific Reports, 12(1), 1-7. [https://doi.org/10.1038/s41598-022-14536-9](https://doi.org/10.1038/s41598-022-14536-9).
* Hershman, R.,* Milshtein, D.,* & Henik, A. (2023). The Contribution of Temporal Analysis of Pupillometry Measurements to Cognitive Research. Psychological Research, 87(1), 28–42. [https://dx.doi.org/10.1007/s00426-022-01656-0](https://dx.doi.org/10.1007/s00426-022-01656-0).
* Hershman, R.,* Share, D. L., Weiss, E. M., Henik. A., & Shechter, A.* (2024). Insights from eye-blinks into the cognitive processes involved in visual word recognition. Journal of Cognition, 7(1): 14, pp. 1–9. [https://doi.org/10.5334/joc.343](https://doi.org/10.5334/joc.343).

### Follow us by adding yourself to our [mailing list](https://cnl.bgu.ac.il/mailing_list/?/register/d5b21996e2446231f719584cfba63b766acad424)

### [https://in.bgu.ac.il/en/Labs/CNL/chap](https://in.bgu.ac.il/en/Labs/CNL/chap)
