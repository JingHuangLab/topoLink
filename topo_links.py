"""
The method to detect topologically linked structure 
of protein complexes from a local perspective.
Code maintained by Yingnan Hou.
"""


from topoly import *
from prody import *
from Bio.PDB import *
import numpy as np
import os
import argparse


def pre_chain_pairs(pdbfile):
    """prepare paired chains of the structure"""
    if pdbfile.endswith("pdb"):
        p = PDBParser()
    else:
        p = MMCIFParser()
    structure = p.get_structure("0", pdbfile)
    chainidlist=[chain.get_id() for chain in structure[0]]
    if " " in chainidlist:
        chainidlist.remove(" ")
    chainpairslist=[]
    # filter out non-complexes:
    if len(chainidlist) < 2:
        return chainpairslist
    for cidi in range(len(chainidlist)-1):
        for cidj in np.arange(cidi+1, len(chainidlist)):
            x=[chainidlist[cidi], chainidlist[cidj]]
            chainpairslist.append(x)
    return chainpairslist


def pre_chain_coords(structure, chainid1, chainid2, discutoff=10,\
    detail=False, outpath='/', ter_rm=15):
    """prepare the coordinates of the paired chains of the structure"""
    chain1=structure.select("chain "+chainid1)
    chain2=structure.select("chain "+chainid2)
    intchain1=structure.select("chain %s and within %f of chain %s"%(chainid1, discutoff, chainid2))
    intchain2=structure.select("chain %s and within %f of chain %s"%(chainid2, discutoff, chainid1))
    
    if (intchain1 == None) or (intchain2 == None):
        return [], [], {}, {}, {}, {}, 0
    ## dict id:resid
    id_c1=chain1.getResindices()
    id_c2=chain2.getResindices()
    resid_c1=chain1.getResnums()
    resid_c2=chain2.getResnums()
    dict_c1={id_c1[i]:resid_c1[i] for i in range(len(id_c1))}
    dict_c2={id_c2[i]:resid_c2[i] for i in range(len(id_c2))}
    
    ## check breaks
    numbreaks=0
    for i in range(min(resid_c1), max(resid_c1)+1):
        if i not in resid_c1:
            numbreaks+=1
    for i in range(min(resid_c2), max(resid_c2)+1):
        if i not in resid_c2:
            numbreaks+=1
    
    ## select interface and combine
    id_intc1=intchain1.getResindices()
    id_intc2=intchain2.getResindices()
    ## check terminals
    ter_rm=max(0,ter_rm)
    min_c1, max_c1=min(id_c1), max(id_c1)
    min_c2, max_c2=min(id_c2), max(id_c2)
    min_id_intc1, max_id_intc1=min(id_intc1), max(id_intc1)
    min_id_intc2, max_id_intc2=min(id_intc2), max(id_intc2)
    if min_id_intc1 - min_c1 <= ter_rm:
        min_id_intc1=min_c1 + ter_rm
    if max_c1 - max_id_intc1 <= ter_rm:
        max_id_intc1=max_c1 - ter_rm
    if min_id_intc2 - min_c2 <= ter_rm:
        min_id_intc2=min_c2 + ter_rm
    if max_c2 - max_id_intc2 <= ter_rm:
        max_id_intc2=max_c2 - ter_rm
    if (max_id_intc1 <= min_id_intc1) or (max_id_intc2 <= min_id_intc2):
        return [], [], {}, {}, {}, {}, 0
    
    selchain1=structure.select("resindex %d to %d"%(min_id_intc1, max_id_intc1))
    selchain2=structure.select("resindex %d to %d"%(min_id_intc2, max_id_intc2))
    selinter=structure.select("resindex %d to %d"%(min_id_intc1, max_id_intc1)+ \
        " or resindex %d to %d"%(min_id_intc2, max_id_intc2))
    id_selc1=selchain1.getResindices()
    id_selc2=selchain2.getResindices()
    num1=len(id_selc1)
    num2=len(id_selc2)
    
    ## loop index
    dict_l1={i+1:id_selc1[i] for i in range(num1)}
    dict_l2={i+1:id_selc2[i] for i in range(num2)}
    c1xyz=selchain1.getCoords().tolist()
    c2xyz=selchain2.getCoords().tolist()
    
    ## if output pdb
    if detail:
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        outfile_pdb=outpath+'/'+"inter_selected_%s%s.pdb"%(chainid1, chainid2)
        prody.writePDB(outfile_pdb, selinter)
    return c1xyz, c2xyz, dict_c1, dict_c2, dict_l1, dict_l2, numbreaks


def identify_pierce_atoms(gln2, dict_c1, dict_l1, threshscore, scanbegin, scanend):
    """identify the piercing atoms"""
    numpoints=len(dict_l1)
    numlink1=0
    maxgln=0
    maxtln=0
    chg1_id_str=''
    minres=round(2*scanend/3)
    if not gln2:
        return numlink1, chg1_id_str, maxgln, maxtln, minres
    elif numpoints <= scanbegin:
        return numlink1, chg1_id_str, maxgln, maxtln, minres
    else:
        scanend= numpoints -1 if scanend >= numpoints else scanend
        gln2matrix=np.array(gln2.get('matrix'))
        max1=np.max(gln2matrix)
        min1=np.min(gln2matrix)
        maxgln= max1 if abs(max1) > abs(min1) else min1
        if abs(maxgln) < (threshscore - 0.01):
            return numlink1, chg1_id_str, maxgln, maxtln, minres
        else:
            chg1_id=[]
            for i in range(numpoints-scanbegin):
                for j in range(scanbegin, scanend):
                    ij1=i+j if i+j < numpoints else numpoints-1
                    maxtln = gln2matrix[i,ij1] if abs(maxtln) < abs(gln2matrix[i,ij1]) else maxtln
                    if abs(gln2matrix[i,ij1]) > threshscore:
                        id1=i+round((ij1-i)/2)
                        resi1=dict_c1.get(dict_l1.get(id1))
                        chg1_id.append(resi1)
                        ## find minimal length to form link
                        lenlink=round((ij1-i+1)/3)
                        minres=lenlink if lenlink < minres else minres

            ## add an opportunity for those almost reach threshold
            if (not chg1_id) and (abs(maxtln) > 0.9*threshscore ):
                for i in range(numpoints-scanend):
                    for j in range(scanend, scanend+6):
                        ij1=i+j if i+j < numpoints else numpoints-1
                        if abs(gln2matrix[i,ij1]) > threshscore:
                            id1=i+round((ij1-i)/2)
                            resi1=dict_c1.get(dict_l1.get(id1))
                            chg1_id.append(resi1)
                            ## find minimal length to form link
                            lenlink=round((ij1-i+1)/3)
                            minres=lenlink if lenlink < minres else minres

            if chg1_id:
                for l in range(min(chg1_id), max(chg1_id)+1):
                    tempstr=''
                    if (l in chg1_id):
                        if l-1 not in chg1_id:
                            tempstr=str(l)+'-'
                            if l-1 not in dict_c1.values():
                                tempstr='--'+tempstr
                        if l+1 not in chg1_id:
                            tempstr=tempstr+str(l)
                            if l+1 in dict_c1.values():
                                numlink1+=1
                                if l < max(chg1_id):
                                    tempstr=tempstr+"|"
                    chg1_id_str=chg1_id_str+tempstr
            return numlink1, chg1_id_str, maxgln, maxtln, minres


def tln_2chains(structure, chainid1, chainid2, discutoff=10, threshscore=0.8,\
    scanbegin=4, scanend=36, detail=False, outpath='./', ter_rm=15):
    """calculate the topological links for the two chains of the structure.
    Return:
        A dictionary with key 'tln' keeps the number of topological links between the
        two chans in the pdbfile, other keys keep information about the whole 
        GLN score, the number of topological links in each chains, the residues index
        of those detected topological links, and the number of missing residues of the
        two chains.
    """
    c1xyz, c2xyz, dict_c1, dict_c2, dict_l1, dict_l2, numbreaks=pre_chain_coords(structure,\
        chainid1, chainid2, discutoff=discutoff, detail=detail, outpath=outpath, ter_rm=ter_rm)
    chainids=chainid1+chainid2
    if not c1xyz:
        result={"tln":0, "wholegln":0, "chain_"+chainid1:0, "chain_"+chainid2:0,\
            "resid_chain_"+chainid1:'', "resid_chain_"+chainid2:'', "res_breaks":numbreaks,\
                "maxgln_"+chainid1:0,"maxgln_"+chainid2:0,\
                    "maxtln_"+chainid1:0,"maxtln_"+chainid2:0, "minres":round(2*scanend/3)}
        return result
    else:
        if detail:
            os.makedirs(outpath, exist_ok=True)
            if len(dict_l1) > 20:
                gln1 = gln(c1xyz, c2xyz, matrix=True,\
                    matrix_map=True,map_filename =outpath+'/'+chainids+"-glnmap-"+chainid2)
                matrix_filename=outpath+'/'+chainids+"-matrix-"+chainid2
                np.save(matrix_filename, np.array(gln1.get('matrix')))
            else:
                gln1 = gln(c1xyz, c2xyz, matrix=True)
            if len(dict_l2) > 20:
                gln2 = gln(c2xyz, c1xyz, matrix=True,\
                    matrix_map=True,map_filename =outpath+'/'+chainids+"-glnmap-"+chainid1)
                matrix_filename=outpath+'/'+chainids+"-matrix-"+chainid1
                np.save(matrix_filename, np.array(gln2.get('matrix')))
            else:
                gln2 = gln(c2xyz, c1xyz, matrix=True)
        else:
            gln1 = gln(c1xyz, c2xyz, matrix=True)
            gln2 = gln(c2xyz, c1xyz, matrix=True)
    numlink1, chg1_id_str, maxgln1, maxtln1, minres1=identify_pierce_atoms(gln2, dict_c1,\
        dict_l1, threshscore, scanbegin, scanend)
    numlink2, chg2_id_str, maxgln2, maxtln2, minres2=identify_pierce_atoms(gln1, dict_c2,\
        dict_l2, threshscore, scanbegin, scanend)
    topolinks=max(numlink1, numlink2)
    minres=min(minres1, minres2)
    wholegln=gln1.get("whole")
    result={"tln":topolinks, "wholegln":wholegln, "chain_"+chainid1:numlink1,\
        "chain_"+chainid2:numlink2, "resid_chain_"+chainid1:chg1_id_str, \
            "resid_chain_"+chainid2:chg2_id_str, "res_breaks":numbreaks,\
                "maxgln_"+chainid1:maxgln1,"maxgln_"+chainid2:maxgln2,\
                    "maxtln_"+chainid1:maxtln1,"maxtln_"+chainid2:maxtln2, "minres":minres}
    return result


def topo_link(pdbfile, outpath='/tmp/topo_links/', detail=False, scanend=36, scanbegin=4,\
    threshscore=0.8, discutoff=10, ter_rm=15):
    """
    Calculates the topological links between chains for the whole structure.
    Parameters:
        pdbfile(str):
            The file of structure to analysis (in *.pdb or *.cif).
        outpath(str, optional):
            The path to save the detailed files.
        detail(bool, optional):
            If the details of the analysis should be output. Default: False.
        scanend(integer, optional):
            The maximum of the scan windows (atoms). Default: 36.
        scanbengin(integer, optional):
            The minimum of the scan windows (atoms). Default: 4.
        threshscore(float, optional):
            The threshold for the absolute value of the GLN scores.Default: 0.8.
        discutoff(float, optional):
            The distance threshhold for the interface selection of the complex.
            Unit is Angstrom. Default: 10.
    Return:
        A nested dictionary with keys of natural numbers (from 1) keep the 
        corresponding values, each of which is the dictionary of the links
        inforfamtion of the paired chains: the key 'tln' keeps the topological
        links number between the paired chains, the key 'wholegln' keeps GLN
        between the two whole chains,other keys keep information about the
        number of topological links in each chains, the residuees index of
        those detected topological links, and the number of missing residues
        of the two chains.
    """
    selename="N CA C"
    structure=prody.parsePDB(pdbfile, subset="bb")
    structure=structure.select("name "+selename)
    chainpairslist=pre_chain_pairs(pdbfile)
    tlnresult={}
    
    if not chainpairslist:
        print("The input data may no be a complex.")
        return tlnresult
    
    for chainpairs in chainpairslist:
        chainid1=chainpairs[0]
        chainid2=chainpairs[1]
        tlnresult[chainid1+chainid2]=tln_2chains(structure, chainid1, chainid2,\
            discutoff=discutoff, threshscore=threshscore, scanbegin=scanbegin,\
                scanend=scanend, detail=detail, outpath=outpath, ter_rm=ter_rm)
    return tlnresult


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', '--pdbfile', type = str, \
        help="The file of structure to analysis (in *.pdb or *.cif).")
    parser.add_argument("-out", "--outpath", type = str,\
        default="/tmp/topo_links/",\
        help="The path to save the output.")
    parser.add_argument("-sf", "--summaryfile", type = str,\
        default="summary_topo_links.txt",\
        help="The file name to save the summaried results of topological links")
    parser.add_argument('--detail', dest='detail', action='store_true',\
        help="If output details")
    parser.add_argument('--no-detail',dest='detail', action='store_false',\
        help="If do not output details")
    parser.set_defaults(detail=False)
    parser.add_argument('-se', "--scanend",type=int, default=36,\
        help="The upper limit of scan windows")
    parser.add_argument('-sb', "--scanbegin",type=int, default=4, \
        help="The lower limit of scan windows")
    parser.add_argument('-ts', "--threshscore",type=float, default=0.8, \
        help="The threshold for the absolute gln scores")
    parser.add_argument('-d', "--discutoff",type=float, default=10.0, \
        help="The distance threshhold for the interface selection of the complex.")
    parser.add_argument('-rm', "--ter_rm",type=int, default=15, \
        help="The length of rediues in terminal to remove for each chain.")
    args = parser.parse_args()

    pdbfile = args.pdbfile
    outpath = args.outpath
    summaryfile=outpath+'/'+args.summaryfile
    detail=args.detail
    scanbegin=args.scanbegin
    scanend=args.scanend
    threshscore=args.threshscore
    discutoff=args.discutoff
    ter_rm=args.ter_rm
    
    os.makedirs(outpath, exist_ok=True)
    if detail:
        filename=os.path.basename( os.path.splitext(pdbfile)[0])
        outpath=outpath+'/'+filename
    x=topo_link(pdbfile, outpath=outpath, detail=detail, scanend=scanend,\
         scanbegin=scanbegin, threshscore=threshscore, discutoff=discutoff, ter_rm=ter_rm)
    os.makedirs(os.path.dirname(summaryfile), exist_ok=True)
    with open(summaryfile, 'a+') as f:
        for chids in x.keys():
            tln=x[chids].get("tln")
            print(pdbfile, chids, tln, x[chids], file=f)


