# Strut-and-Tie-automated-modeling
The codes for the STM modeling investigations

## Demo case based on the paper: 
- Reference: Xia, Y., Langelaar, M., & Hendriks, M. A. (2020). A critical evaluation of topology optimization results for strut‐and‐tie modeling of reinforced concrete. Computer‐Aided Civil and Infrastructure Engineering, 35(8), 850-869.
  
- Run Case1_beam.m in the file of "Extraction and analyzing demo".
  
- Contents:
1. TO result extraction method for transforming the optimized topology into a truss-like structure.
2. Evaluation of topology optimization results for strut-and-tie modeling of reinforced concrete.
3. Input: Topology optimization results.
3. Output: Truss-like structure and three evaluation indexes. TRS,SR and STS index are three indexes for evaluating truss-like models.
    - TRS index: Tensile region similarity index, tensile stress fields of the original structure and TO results are compared.
    - SR index: Steel reinforcement ratio, the SR ratio is calculated as the volume fraction of steel with respect to the concrete volume.
    - STS index: Suitable truss structure index, it measures the degree to which the obtained truss-like structure.
