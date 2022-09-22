# Topological links detection

The method to efficiently detect topologically linked structure of protein-protein complexes from a local perspective.

This program is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Requirements

* Python 3 (3.5 or later)

#### Setup Conda Environment (Optional)

1. Create a new conda environment (eg: topolinks):  ```conda create -n topolinks python=3.7 ```
2. Activate this environment: ``` conda activate topolinks```

#### Install dependent packages

*   **topoly**: pip install topoly
*   **ProDy**: pip install ProDy

## Running

*  **For one structure**
  
1. Run ```topo_links.py``` via ```-in``` parameter pointing to the file of protein complex (in *.pdb or *.cif)  you wish to analysis. The simplest way to test it with an example:
   ```python topo_links.py  -in ./examples/glnB-gspJ.pdb```

2. More choices to run ```topo_links.py``` through the follow parameters:
   * -out: (str, optional) the output directory to save the outputs. Defalut: ```/tmp/topo_links```
   * -sf: (str, optional) the file name of the summary of the results. Defalut: ```summary_topo_links.txt```
   * --detail: if output the details of the analysis.
   * --no-detail: if do not output the details of the analysis. (Default)
   * -se: (int, optional) the upper limit of scan windows. Default: 36.
   * -sb: (int, optional) the lower limit of scan windows. Default: 4.
   * -ts: (float, optional) the threshold score for the alsulate gln scores. Default: 0.8.
   * -d: (float, optional) the distance threshold for the selection of interface of the complexes. Default: 10.0.
   * -rm: (int, optional) the length of residues in terminal to remove for each chain. Default: 15.

*  **Batch analysis for multiple structures**

1. Run ```batch_topo_links.sh``` pointing to the folder where the structures of protein complexes are located. Note that all the structures in  *.pdb or *.cif in the folder will be analyzed. The simplest way to test it with an example: 
  ```./batch_topo_links.sh ./examples/```

2. More choices to run ```batch_topo_links.sh``` through the arguments which are refered to the positional parameters.
   
   * ```$1``` (the first argument) refers to the  the folder where the structures of protein complexes are located.
   * ```$2``` (the second argument, optional) refers to the output directory to save the output. Defalut: ```/tmp/topo_links```.
   * ```$3``` (the third argument, optional) refers to the file name of the summary of the results. Defalut: ```summary_topo_links.txt```
   * ```$4``` (the fourth argument, optional) refers to whether output the details of the analysis. Use ```--detail``` to  output the details of the analysis. And ```--no-detail``` refer to do not output the details of the analysis (Default).
   *  ```$5``` (the fifth argument, optional) refers to the upper limit of scan windows. Default: 36.
   *  ```$6``` (the sixth argument, optional) refers to the lower limit of scan windows. Default: 4.
   *  ```$7``` (the seventh argument, optional) refers to the threshold score for the alsulate gln scores. Default: 0.8.
   *  ```$8``` (the eighth argument, optional) refers to the distance threshold for the selection of interface of the complexes. Default: 10.0.
   *  ```$9``` (the ninth argument, optional) refers to the length of residues in terminal to remove for each chain. Default: 15.


## Output

* The outputs will be saved in the directory provided via the ```-out``` parameter of the ```topo_links.py```, or the second argument of the ```batch_topo_links.sh```. Defalut: ```/tmp/topo_links```.
* The outputs include the file of summaried results named via ```-sf``` in ```topo_links.py``` or  the third argument of the ```batch_topo_links.sh``` (Defalut: "summary_topo_links.txt"). Four columns are in the summaried file: the first column is the file of the structure, the second cloumn is the paired chains of the structure, the third column is the number of topological links between the two paired chains and the fourth column is the nested dictionary strored the information of topological links between the two paired chains.
* If ```--detail``` in running ```topo_links.py``` or the fourth argument in runing ```batch_topo_links.sh```, a folder contains the selected interface of paired chains, the GLN matrix and the corresponding plot for each of the paired chain.



## Reference

Yingnan Hou, Tengyu Xie, Liuqing He, Liang Tao, Jing Huang. Topological Links in Predicted Protein Complex Structures Reveal Limitations of AlphaFold. bioRxiv 2022.09.14.507881; doi: https://doi.org/10.1101/2022.09.14.507881.



## Acknowledgments

* We would like to thank the [Topoly](https://topoly.cent.uw.edu.pl/#) team for developing the topoly package.
* We would like to thank the [ProDy](http://prody.csb.pitt.edu/index.html) team for developing the ProDy package.
* We would like to thank the [AlphaFold](https://github.com/deepmind/alphafold) team for the excellent job on protein structure prediction.

