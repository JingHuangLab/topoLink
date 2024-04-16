# Topological links detection

The method to efficiently detect topologically linked structure of protein-protein complexes from a local perspective.

This program is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Requirements

* Python 3 (3.5 or later)
* topoly (between 0.9.24 and 1.0.0)
* ProDy


#### Setup Conda Environment (Optional)

1. Create a new conda environment (eg: topolinks):  ```conda create -n topolinks python=3.7 ```
2. Activate this environment: ``` conda activate topolinks```

#### Install dependent packages

*   **topoly**: pip install topoly==1.0.0
*   **ProDy**: pip install ProDy


## How to run

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


## Interpreting the results

#### summary

* For using ```topo_links.py```, the path and file name of the summaried results are specified by the ```-out``` and ```-sf``` parameters respectively. For using ```batch_topo_links.sh```, the path and file name of the summaried results are specified by the second and third arguments respectively. By default, you can find the results at ```/tmp/topo_links/summary_topo_links.txt```.
  
Let's take the summaried results of ```./examples/glnB-gspJ.pdb``` as an example. Four columns are in the summary file as following: 
```./examples/glnB-gspJ.pdb AB 1 {'tln': 1, 'wholegln': -0.943, 'chain_A': 1, 'chain_B': 1, 'resid_chain_A': '48-54', 'resid_chain_B': '94-98', 'res_breaks': 0, 'maxgln_A': -1.0, 'maxgln_B': -0.986, 'maxtln_A': -0.906, 'maxtln_B': -0.978, 'minres': 6}```

1. ```./examples/glnB-gspJ.pdb``` is the first column (the file of the structure).
2. ```AB``` is the second column (the paired chains of the structure, i.e., chain A and chain B).
3. ```1``` is the third column (the number of topological links between chain A and chain B). **By checking whether the value in the third column is greater than 0, you can quickly determine whether there is a topological link between the two corresponding chains.** Here we can see that there is a topological link between chain A and chain B.
4. Others is the fourth column (Detailed information stored in python dictionary format).
   * ```tln``` is the **t**opological **l**ink **n**umber between chain A and chain B (here equals to 1).
   * ```wholegln``` is the [GLN](https://en.wikipedia.org/wiki/Linking_number) value of the whole chains A and B at the interface (here equals to -0.943).
   * ```chain_A``` is the number of topological links detected in chain A (here equals to 1).
   * ```chain_B``` is the number of topological links detected in chain B (here equals to 1).
   * ```resid_chain_A``` is the residue indexes where the detected topological links in chain A. ```48-54``` means that the fragement between residues 48 to 54 in chain A contributes to the topological link. When there are multiple discontinuous fragments contributes to the topological links, they are separated by "|".
   * ```resid_chain_B``` is the residue indexes where the detected topological links in chain B. ```94-98``` means that the fragement between residues 94 to 98 in chain B contributes to the topological link.
   * ```res_breaks``` is the number of missing residues in chain A and chain B at the interface (Note that the broken chains may cause fake topological links).
   * ```maxgln_A``` is the GLN value with the largest absolute value on chain A (see ```AB-matrix-A.npy``` if ```--detail``` is enabled), where the GLN value is calculated for all possible fragments on chain A with the whole chain B.
   * ```maxgln_B``` is the GLN value with the largest absolute value on chain B (see ```AB-matrix-B.npy``` if ```--detail``` is enabled), where the GLN value is calculated for all possible fragments on chain B with the whole chain A.
   * ```maxtln_A``` is the GLN value with the largest absolute value on chain A, where the GLN values are calculated for the fragments _whose length within ```-sb``` to ```-se```_ on chain A with the whole chain B.
   * ```maxtln_B``` is the GLN value with the largest absolute value on chain B, where the GLN values are calculated for the fragments _whose length within ```-sb``` to ```-se```_ on chain B with the whole chain A.
   * When a topological link is detected, ```minres``` is the length (in number of residues) of the shortest fragment that is judged to form the topological link. When there is no topological link detected, it is the default value ($2*se/3$).



#### Detail

* If ```--detail``` in running ```topo_links.py``` or the fourth argument in runing ```batch_topo_links.sh```, a folder named by the analyzed file name will be generated, containing files as following (e.g. a protein contains chain A and chain B):
  
  * The file ```inter_selected_AB.pdb``` stores the isolated interface for chain A and chain B. All chain A and chain B (including the whole chain A and the whole chain B) in this README refer to their parts at the interface.
  * For ```AB-glnmap-A.png```, ```AB``` in the file name means that the paired chains are chain A and chain B; ```-A``` in the file name means that this glnmap is for chain A. For this figure, each point of the heatmap represents a fragment of chain A (the row represents the starting atom index, and the column represents the ending atom index). The color represents the GLN value between this fragment and chain B. The larger the absolute GLN value, the more intense the degree of entanglement. By default, when the absolute value is greater than 0.8, the corresponding fragment will be judged to have a topological link. Similar, ```AB-glnmap-B.png``` means that the heatmap of the GLN values for chain B.
  * The ```AB-matrix-A.npy``` file stores the GLN value corresponding to the figure ```AB-glnmap-A.png```. The ```AB-matrix-B.npy``` file stores the GLN value corresponding to the figure ```AB-glnmap-B.png```.


Note that if the protein-protein complex contains more than two chains, all pairwise combinations of these chains will be analyzed.


## Reference

Yingnan Hou, Tengyu Xie, Liuqing He, Liang Tao, and Jing Huang*. Topological Links in Predicted Protein Complex Structures Reveal Limitations of AlphaFold. *Communications Biology*, 2023, 6(1): 1098. DOI: https://doi.org/10.1038/s42003-023-05489-4
  



## Acknowledgments

* We would like to thank the [Topoly](https://topoly.cent.uw.edu.pl/#) team for developing the topoly package.
* We would like to thank the [ProDy](http://prody.csb.pitt.edu/index.html) team for developing the ProDy package.
* We would like to thank the [AlphaFold](https://github.com/deepmind/alphafold) team for the excellent job on protein structure prediction.

