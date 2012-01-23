import models.sxrd_test5_sym_new_test_new66_2 as model
from models.utils import UserVars
import numpy as np
from operator import mul
from numpy.linalg import inv
import sys
sys.path.append('Y:\\codes\\my code\\modeling files\\surface modeling 1\\scripts')
import domain_creator2

####################################################################
###########to be changed in the grid settings if domains have changed###############
#1 missing Fe layers, maybe not appicable, maybe missed at a different layer
#2 set the reference point for extra lone oxygen pair, z value ranges should be reset
#3 setFeu_n will be changed possibly, delete or add according to the situation
unitcell = model.UnitCell(5.038, 5.434, 7.3707, 90, 90, 90)
inst = model.Instrument(wavel = .833, alpha = 2.0)
bulk = model.Slab(T_factor='B')
domain0 =  model.Slab(c = 1.0,T_factor='B')

bulk.add_atom( "Fe2", "Fe", 0.00000e+00 ,     8.30000e-01 ,     8.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe3", "Fe", 5.00000e-01 ,     3.30000e-01 ,     8.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe4", "Fe", 5.00000e-01 ,     8.80000e-01 ,     6.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe6", "Fe", 0.00000e+00 ,     3.79000e-01 ,     6.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe8", "Fe", 0.00000e+00 ,     7.61000e-01 ,     3.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe9", "Fe", 5.00000e-01 ,     2.60000e-01 ,     3.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe10", "Fe", 5.00000e-01 ,     8.10000e-01 ,     1.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "Fe12", "Fe", 0.00000e+00 ,     3.10000e-01 ,     1.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O1", "O",  6.53000e-01 ,     9.73000e-01 ,     9.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O2", "O",  8.47000e-01 ,     4.73000e-01 ,     9.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O3", "O",  3.06000e-01 ,     6.05000e-01 ,     7.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O4", "O",  1.94000e-01 ,     1.04000e-01 ,     7.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O5", "O",  8.47000e-01 ,     7.37000e-01 ,     5.97000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O6", "O",  6.53000e-01 ,     2.36000e-01 ,     5.97000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O7", "O",  3.47000e-01 ,     9.04000e-01 ,     4.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O8", "O",  1.53000e-01 ,     4.03000e-01 ,     4.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O9", "O",  6.94000e-01 ,     5.35000e-01 ,     2.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O10", "O",  8.06000e-01 ,     3.50000e-02 ,     2.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O11", "O",  1.53000e-01 ,     6.67000e-01 ,     9.70000e-02 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
bulk.add_atom( "O12", "O",  3.47000e-01 ,     1.67000e-01 ,     9.70000e-02 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )
#domain0 here is a reference domain, the atoms are ordered according to hight (z values)
#it is a super surface structure by stacking the surface slab on bulk slab, the repeat vector was counted 
domain0.add_atom( "O1_1_0", "O",  6.53000e-01 ,     1.11210e+00 ,     1.90300e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_2_0", "O",  8.47000e-01 ,     6.12100e-01 ,     1.90300e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_2_0", "Fe", 0.00000e+00 ,     9.69100e-01 ,     1.85500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1. )    
domain0.add_atom( "Fe1_3_0", "Fe", 5.00000e-01 ,     4.69100e-01 ,     1.85500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1. )    
domain0.add_atom( "O1_3_0", "O",  3.06000e-01 ,     7.44100e-01 ,     1.75000e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_4_0", "O",  1.94000e-01 ,     2.43100e-01 ,     1.75000e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_4_0", "Fe", 5.00000e-01 ,     1.01910e+00 ,     1.64500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_6_0", "Fe", 0.00000e+00 ,     5.18100e-01 ,     1.64500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_5_0", "O",  8.47000e-01 ,     8.76100e-01 ,     1.59700e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_6_0", "O",  6.53000e-01 ,     3.75100e-01 ,     1.59700e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_7_0", "O",  3.47000e-01 ,     1.04310e+00 ,     1.40300e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_8_0", "O",  1.53000e-01 ,     5.42100e-01 ,     1.40300e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_8_0", "Fe", 0.00000e+00 ,     9.00100e-01 ,     1.35500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_9_0", "Fe", 5.00000e-01 ,     3.99100e-01 ,     1.35500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_9_0", "O",  6.94000e-01 ,     6.74100e-01 ,     1.25000e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_10_0", "O",  8.06000e-01 ,     1.74100e-01 ,     1.25000e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_10_0", "Fe", 5.00000e-01 ,     9.49100e-01 ,     1.14500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe1_12_0", "Fe", 0.00000e+00 ,     4.49100e-01 ,     1.14500e+00 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_11_0", "O",  1.53000e-01 ,     8.06100e-01 ,     1.09700e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O1_12_0", "O",  3.47000e-01 ,     3.06100e-01 ,     1.09700e+00 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    

domain0.add_atom( "O1_0", "O",  6.53000e-01 ,     9.73000e-01 ,     9.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O2_0", "O",  8.47000e-01 ,     4.73000e-01 ,     9.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe2_0", "Fe", 0.00000e+00 ,     8.30000e-01 ,     8.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1. )    
domain0.add_atom( "Fe3_0", "Fe", 5.00000e-01 ,     3.30000e-01 ,     8.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1. )    
domain0.add_atom( "O3_0", "O",  3.06000e-01 ,     6.05000e-01 ,     7.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O4_0", "O",  1.94000e-01 ,     1.04000e-01 ,     7.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe4_0", "Fe", 5.00000e-01 ,     8.80000e-01 ,     6.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe6_0", "Fe", 0.00000e+00 ,     3.79000e-01 ,     6.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O5_0", "O",  8.47000e-01 ,     7.37000e-01 ,     5.97000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O6_0", "O",  6.53000e-01 ,     2.36000e-01 ,     5.97000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O7_0", "O",  3.47000e-01 ,     9.04000e-01 ,     4.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O8_0", "O",  1.53000e-01 ,     4.03000e-01 ,     4.03000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe8_0", "Fe", 0.00000e+00 ,     7.61000e-01 ,     3.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe9_0", "Fe", 5.00000e-01 ,     2.60000e-01 ,     3.55000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O9_0", "O",  6.94000e-01 ,     5.35000e-01 ,     2.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O10_0", "O",  8.06000e-01 ,     3.50000e-02 ,     2.50000e-01 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe10_0", "Fe", 5.00000e-01 ,     8.10000e-01 ,     1.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "Fe12_0", "Fe", 0.00000e+00 ,     3.10000e-01 ,     1.45000e-01 ,     3.20000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O11_0", "O",  1.53000e-01 ,     6.67000e-01 ,     9.70000e-02 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    
domain0.add_atom( "O12_0", "O",  3.47000e-01 ,     1.67000e-01 ,     9.70000e-02 ,     3.30000e-01 ,     1.00000e+00 ,     1.00000e+00 )    

#id list according to the order in the reference domain
ref_id_list=["O1_1_0","O1_2_0","Fe1_2_0","Fe1_3_0","O1_3_0","O1_4_0","Fe1_4_0","Fe1_6_0","O1_5_0","O1_6_0","O1_7_0","O1_8_0","Fe1_8_0","Fe1_9_0","O1_9_0","O1_10_0","Fe1_10_0","Fe1_12_0","O1_11_0","O1_12_0",\
"O_1_0","O_2_0","Fe_2_0","Fe_3_0","O_3_0","O_4_0","Fe_4_0","Fe_6_0","O_5_0","O_6_0","O_7_0","O_8_0","Fe_8_0","Fe_9_0","O_9_0","O_10_0","Fe_10_0","Fe_12_0","O_11_0","O_12_0"]
#the matching row Id information in the symfile
sym_file_Fe=np.array(['Fe1_0','Fe2_0','Fe3_0','Fe4_0','Fe5_0','Fe6_0','Fe7_0','Fe8_0','Fe9_0','Fe10_0','Fe11_0','Fe12_0',\
    'Fe1_1_0','Fe1_2_0','Fe1_3_0','Fe1_4_0','Fe1_5_0','Fe1_6_0','Fe1_7_0','Fe1_8_0','Fe1_9_0','Fe1_10_0','Fe1_11_0','Fe1_12_0'])
sym_file_O=np.array(['O1_0','O2_0','O3_0','O4_0','O5_0','O6_0','O7_0','O8_0','O9_0','O10_0','O11_0','O12_0',\
    'O1_1_0','O1_2_0','O1_3_0','O1_4_0','O1_5_0','O1_6_0','O1_7_0','O1_8_0','O1_9_0','O1_10_0','O1_11_0','O1_12_0'])
batch_path_head='Y:\\codes\\my code\\modeling files\\new domain test\\genx file\\good ones\\scripts\\batch_files\\'
sym_file_head='Y:\\codes\\my code\\modeling files\\surface modeling 1\\scripts\\'
#create a domain class and initiate the chemical equivalent domains
#when change or create a new domain, make sure the terminated_layer set right
#domain1
rgh_domain1=UserVars()
domain_class_1=domain_creator2.domain_creator(ref_domain=domain0,id_list=ref_id_list,terminated_layer=0,domain_N=1,new_var_module=rgh_domain1)
domain1A=domain_class_1.domain_A
domain1B=domain_class_1.domain_B
#domain2
rgh_domain2=UserVars()
domain_class_2=domain_creator2.domain_creator(ref_domain=domain0,id_list=ref_id_list,terminated_layer=4,domain_N=2,new_var_module=rgh_domain2)
domain2A=domain_class_2.domain_A
domain2B=domain_class_2.domain_B
#add sorbates for two domains
#when change the domain, only update the attach_atm_id info
domain_class_1.add_sorbate_polyhedra(domain=domain1A,polyhedra_flag='hexahedra',extra_flag='1_1+0_1',extra_flag2='type1',attach_atm_id=[['O1_1_0','O1_2_0'],['O1_3_0','O1_4_0']],el='Pb',id_attach=['A','AA'])
domain_class_1.add_sorbate_polyhedra(domain=domain1B,polyhedra_flag='hexahedra',extra_flag='1_1+0_1',extra_flag2='type1',attach_atm_id=[['O1_7_0','O1_8_0'],['O1_9_0','O1_10_0']],el='Pb',id_attach=['B','BB'])
domain_class_2.add_sorbate_polyhedra(domain=domain2A,polyhedra_flag='hexahedra',extra_flag='1_1+0_1',extra_flag2='type1',attach_atm_id=[['O1_5_0','O1_6_0'],['O1_7_0','O1_8_0']],el='Pb',id_attach=['A','AA'])
domain_class_2.add_sorbate_polyhedra(domain=domain2B,polyhedra_flag='hexahedra',extra_flag='1_1+0_1',extra_flag2='type1',attach_atm_id=[['O1_11_0','O1_12_0'],['O1_0','O2_0']],el='Pb',id_attach=['B','BB'])

#set domain1A domain1B to class1
domain_class_1.domain1A=domain1A
domain_class_1.domain1B=domain1B
domain_class_2.domain2A=domain2A
domain_class_2.domain2B=domain2B
#set new variables
#be carefule about the N_list, the first two items will always change for different terminated surface domain
domain_class_1.set_new_vars(head_list=['u_o_n','u_Fe_n','dx_n','dy_n','dz_n','oc_n','dx_sign_n','dy_sign_n','dz_sign_n'],N_list=[4,3,7,7,7,7,7,7,7])
domain_class_2.set_new_vars(head_list=['u_o_n','u_Fe_n','dx_n','dy_n','dz_n','oc_n','dx_sign_n','dy_sign_n','dz_sign_n'],N_list=[5,2,7,7,7,7,7,7,7])
#some other parameters to be used 
#should be the same
domain_class_1.set_discrete_new_vars_batch(batch_path_head+'new_varial_file_hexahedra_class1.txt')
domain_class_2.set_discrete_new_vars_batch(batch_path_head+'new_varial_file_hexahedra_class2.txt')
#do grouping for top seven layers
#when change domain,just update the first_atom_id info
atm_gp_list_domain1=domain_class_1.grouping_sequence_layer(domain=[domain1A,domain1B], first_atom_id=['O1_1_0','O1_7_0'],\
    sym_file={'Fe':sym_file_head+'Fe0 output file for Genx reading.txt','O':sym_file_head+'O0 output file for Genx reading.txt'},\
    id_match_in_sym={'Fe':sym_file_Fe,'O':sym_file_O},layers_N=7,use_sym=True)
atm_gp_list_domain2=domain_class_2.grouping_sequence_layer(domain=[domain2A,domain2B], first_atom_id=['O1_5_0','O1_11_0'],\
    sym_file={'Fe':sym_file_head+'Fe0 output file for Genx reading.txt','O':sym_file_head+'O0 output file for Genx reading.txt'},\
    id_match_in_sym={'Fe':sym_file_Fe,'O':sym_file_O},layers_N=7,use_sym=True)
#SET the atom group list to class varials,note the attribute can have the same name in different class
domain_class_1.atm_gp_list_domain1=atm_gp_list_domain1
domain_class_2.atm_gp_list_domain1=atm_gp_list_domain2
#the first atom group will be the reference group for scaling operation of dx dy dz
ref_atm_gp_domain1=atm_gp_list_domain1[0]
ref_atm_gp_domain2=atm_gp_list_domain2[0]
#set to class variables, still note the same name for different class
domain_class_1.ref_atm_gp_domain1=ref_atm_gp_domain1
domain_class_2.ref_atm_gp_domain1=ref_atm_gp_domain2
#group the sorbate of Pb and oxygen pair
#txt file kept the same
atm_gp_Pb_domain1,atm_gp_Pb2_domain1,atm_gp_Os_1_domain1,atm_gp_Os_2_domain1=domain_class_1.grouping_discrete_layer_batch(batch_path_head+'group_discrete_layer_file_hexahedra_class1.txt')
atm_gp_Pb_domain2,atm_gp_Pb2_domain2,atm_gp_Os_1_domain2,atm_gp_Os_2_domain2=domain_class_2.grouping_discrete_layer_batch(batch_path_head+'group_discrete_layer_file_hexahedra_class2.txt')
#make a domain libratry wrapping two chemical equivalent domains
domain={'domain1A':{'slab':domain1A,'wt':1.},'domain1B':{'slab':domain1B,'wt':0.},'domain2A':{'slab':domain2A,'wt':0.},'domain2B':{'slab':domain2B,'wt':0.}}
sample = model.Sample(inst, bulk, domain, unitcell,coherence=False,surface_parms={'delta1':0.,'delta2':0.1391})
#move on
def Sim(data):
#scale class one batch
#extract the fitting par values in the associated attribute and then do the scaling(initiate+processing)
#in sim_batch file, ocu lines should be changed accordingly, and make sure the symbol used in txt is the same as that defined in the class
#note the atomlist and ref_value has the same name even at the different class
#in scale txt file, change the index order accordingly
    domain_class_1.init_sim_batch(batch_path_head+'sim_batch_class1.txt')
    domain_class_1.scale_opt_batch(batch_path_head+'scale_operation_file_class1.txt')
#scale class two
    domain_class_2.init_sim_batch(batch_path_head+'sim_batch_class2.txt')
    domain_class_2.scale_opt_batch(batch_path_head+'scale_operation_file_class2_1.txt')
#change the orientation of the initial polyhedra
    domain_class_1.updata_polyhedra_orientation_batch(batch_path_head+'hexahedra_orientation_batch_file_class1.txt')
    domain_class_2.updata_polyhedra_orientation_batch(batch_path_head+'hexahedra_orientation_batch_file_class2.txt')    
#updata the body center point of each polyhedra
    domain_class_1.updata_polyhedra_center_point_batch(batch_path_head+'hexahedra_center_point_batch_file_class1.txt')
    domain_class_2.updata_polyhedra_center_point_batch(batch_path_head+'hexahedra_center_point_batch_file_class2.txt')
#updata sorbate xyz (bidentate configuration here)
#this part should kept the same
    domain_class_1.updata_sorbate_polyhedra2_batch(batch_path_head+'update_sorbate_hexahedra_batch_file_class1.txt')
    domain_class_2.updata_sorbate_polyhedra2_batch(batch_path_head+'update_sorbate_hexahedra_batch_file_class2.txt')
#updata and normalize domain weight values
    wt_tot=rgh_domain1.domain_wt1+rgh_domain1.domain_wt2+rgh_domain2.domain_wt3+rgh_domain2.domain_wt4
    domain['domain1A']['wt']=rgh_domain1.domain_wt1/wt_tot
    domain['domain1B']['wt']=rgh_domain1.domain_wt2/wt_tot
    domain['domain2A']['wt']=rgh_domain2.domain_wt3/wt_tot
    domain['domain2B']['wt']=rgh_domain2.domain_wt4/wt_tot
#9.a loop through the data sets
    F =[]
    for data_set in data:
        # 9.b create all the h,k,l values for the rod (data_set)
        h = data_set.extra_data['h']
        k = data_set.extra_data['k']
        l = data_set.x
        # 9.c. calculate roughness using beta model
        #LB = data_set.extra_data['LB']
        #dL = data_set.extra_data['dL']
        #rough = (1-beta)/((1-beta)**2 + 4*beta*np.sin(np.pi*(l - LB)/dL)**2)**0.5
        # 9.d. Calculate the structure factor
        f = sample.calc_f(h, k, l)
        # 9.e Calculate |F|
        i = abs(f)
        # 9.f Append the calculated intensity to the list I
        F.append(i)
    return F