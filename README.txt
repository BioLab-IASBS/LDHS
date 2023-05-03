 Harmony search algorithm for detecting interactions with grouped and weighted SNPs
 
This study proposes a new approach in order to improve swarm intelligence based algorithms such as harmony search for finding disease-causing SNP interactions. 
The idea is applied on MP-HS-DHSI algorithm, and therefore the available code for this research is used (and modified) in order to evaluate the result. 
======================================================================

The flowchart of the proposed method is represented in Figure 1 of the paper. 
(Flowchart of steps taken in this research with a focus on initializing harmony memories)

  

  -----
>> we proposed a SNP grouping method  based on the introduced SNP categories in this section, i.e. the SNPs of the coding regions, the SNPs of the promoter regions, and the SNPs of the non-coding regions, and also considering the location of SNPs. To be more exact, this grouping provides sets of SNPs that are located in DNA regions that are naturally related to one goal – the expression of a specific gene and finally its function as a protein. The SNP information (location and category) is available in the SNP database on NCBI website .  
In addition to considering different types of SNPs (based on their location)  in developing the new idea, an important biological knowledge is that genetic contents that are located close to each other on DNA sequence (on the same chromosome) have more tendency to be inherited together.  In conclusion  , SNPs located on genes and promoters can be more important (disease-causing) that SNPs located on introns.   
These SNPs are located close to each other, and therefore have a high chance of being inherited together. 
After retrieving the SNP information and grouping them, the number of vectors (harmonies or individual solutions) generated using each group is determined according to the number of functional or non-functional SNPs in each group. The number of individuals that are generated using functional SNPs is more (experimentally specified) since functional SNPs have more chance to cause disease. In other words, a weight is assigned to each group of SNPs. This weight determines in generating how many individuals in the HM this group should be considered. Algorithm initial population SNPs are generated using grouped and weighted SNPs in the file HS_2019_multiCRITERIA5.


 
code list:

(A modified version of MP-HS-DHSI)
  ============
     |-- NEW_DME.m  : the main program for DME models  (improved version of MP_HS_DHSI_for_DME.m*)
    
     |-- NEW_NDME.m : the main program for NDME mdoels  (improved version of MP_HS_DHSI_for_NDME.m*)
   
     |-- NEW_HS.m : the search program with weighted and grouped SNPs. ((improved version of HS_2019_multiCRITERIA5.m*)
   
     |-- multi_criteria 3.m : the score functions for evaluating the association (no change compared to MP-HS-DHSI* original implementation)
   
     |-- Gtest_score.m : the G-test function (no change compared to original code from MP-HS-DHSI*)
   
     |-- MDR2.m  : The MDR program for the 3rd stage (no change compared to original code from MP-HS-DHSI*)


* Tuo, S., Liu, H. and Chen, H., 2020. Multipopulation harmony search algorithm for the detection of high-order SNP interactions. Bioinformatics, 36(16), pp.4389-4398.
Code related to MP-HS-DHSI is downloaded from https://github.com/shouhengtuo/MP-HS-DHSI.


