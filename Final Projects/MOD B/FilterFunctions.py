import ROOT
import ROOT.ROOT as rr

import cpp

dirOutPath = '/data/Skim/'

# - Dictionary for iSkim variables and labels #####################################

default_nbins = 30

SkimRanges = {
    2 : 
                {'Muon_pt15':(default_nbins, 0, 200),
                 'mu_pt': (default_nbins, 0, 200),
                 'el_pt': (default_nbins, 0, 200),
                 'max_eta': (default_nbins, 0, 4),
                 'inv_m': (default_nbins, 0, 200),
                 'MET_pt': (default_nbins, 0, 200),
                 'dRHJ': (default_nbins, 0, 4),
                 'dRLJ': (default_nbins, 0, 4),
                 'ST': (default_nbins, 0, 2000),
                 'dPhi0': (default_nbins, -4, 4),
                 'dPhi1': (default_nbins, -4, 4),},
    
    3 :
                {'Muon_pt15':(default_nbins, 0, 200),
                 'mu_pt0': (default_nbins, 0, 200),
                 'inv_m01': (default_nbins, 0, 200),
                 'inv_m12': (default_nbins, 0, 200),
                 'inv_m02': (default_nbins, 0, 200),   
                 'inv_m3': (default_nbins, 0, 200),
                 'ST': (default_nbins, 200, 2000),
                 'Jet_pt0': (default_nbins, 0, 200),
                 'bJet_pt0': (default_nbins, 0, 200),
                 'MET_pt': (default_nbins, 0, 200),
                 'dPhi0': (default_nbins, -4, 4),
                 'dPhi1': (default_nbins, -4, 4),
                 'dPhi2': (default_nbins, -4, 4),
                 'dR01' : (default_nbins, 0, 4),
                 'dR02' : (default_nbins, 0, 4),}
                 #'di_Jet_invm': (default_nbins, 0, 200),}       
        }

SkimLabels = { 
    2 : 
                {'Muon_pt15': 'Muon p_{T} 15 / GeV',
                 'mu_pt': 'Muon p_{T} / GeV',
                 'el_pt': 'Electron p_{T} / GeV',
                 'max_eta': 'Max #eta',
                 'inv_m': 'Inv m / GeV',
                 'MET_pt': 'Missing E',
                 'dRHJ': 'dRHJ',
                 'dRLJ': 'dRLJ',
                 'ST': 'ST',
                 'dPhi0': 'dPhi0',
                 'dPhi1': 'dPhi1',},
    
    3 :
                {'Muon_pt15': 'Muon p_{T} 15 / GeV',
                'mu_pt0': 'Muon p_{T} / GeV',
                'inv_m01': 'Inv m 01',
                'inv_m12': 'Inv m 12',
                'inv_m02': 'Inv m 02',
                'inv_m3': 'Inv m 3',
                'ST': 'ST',
                'Jet_pt0': 'Jet pt0',
                'bJet_pt0': 'b-Jet pt0',
                'MET_pt': 'Missing E',
                'dPhi0': 'dPhi0',   
                'dPhi1': 'dPhi1',   
                'dPhi2': 'dPhi2',
                'dR01' : 'dR01',
                'dR02' : 'dR02',}
                #'di_Jet_invm': 'di_Jet_invm',}
        }

# - iSkim1 ########################################################################
def FSkim1(df):
    fdf = df.Filter('iSkim == 1', 'iSkim1')\
            .Define('maskMu', 'Muon_pfRelIso03_all < 0.15 && Muon_pt > 7')\
            .Define('nGoodMu', 'Sum(maskMu)')\
            .Define('Muon_pt15', 'Muon_pt[maskMu]')\
            .Define('maskEl', 'Electron_pfRelIso03_all < 0.15')\
            .Define('nGoodEl', 'Sum(maskEl)')\
            .Filter('nGoodMu == 2 && nGoodEl == 0', 'GoodEvent')\
            .Define('maskJClean', 'cleanJ(Jet_eta, Muon_eta[maskMu], Jet_phi, Muon_phi[maskMu])')\
            .Define('maskBjet', 'Jet_btagDeepB > 0.7')\
            .Define('nJClean', 'Sum(maskJClean)')\
            .Define('maskBJClean', 'prod(maskJClean,maskBjet)')\
            .Define('nBJ', 'Sum(maskBJClean)')\
            .Filter('nJClean >= 4 && nBJ > 0 && nBJ < 2', 'GoodJet') 
            
    return fdf
# .Filter('nCleanJet >= 4 && nBJet > 0 && nBJet <= 2', 'GoodJet')

# - iSkim2 ########################################################################
def FSkim2(df):
    fdf = df.Filter('iSkim == 2', 'iSkim2')\
            .Define('maskMu', 'Muon_pfRelIso03_all < 0.1 && Muon_pt > 27')\
            .Define('nGoodMu', 'Sum(maskMu)')\
            .Define('Muon_pt15', 'Muon_pt[maskMu]')\
            .Define('maskEl', 'Electron_pfRelIso03_all < 0.1 && Electron_pt > 20')\
            .Define('nGoodEl', 'Sum(maskEl)')\
            .Define('Electron_pt15', 'Electron_pt[maskEl]')\
            .Filter('nGoodMu == 1 && nGoodEl == 1', 'GoodEvent')\
            .Define('maskJClean', 'cleanJ(Jet_eta, Merge(Muon_eta[maskMu], Electron_eta[maskEl]),\
                    Jet_phi,Merge(Muon_phi[maskMu], Electron_phi[maskEl]))')\
            .Define('maskBjet', 'Jet_btagDeepB > 0.85')\
            .Define('nJClean', 'Sum(maskJClean)')\
            .Define('maskBJClean', 'prod(maskJClean,maskBjet)')\
            .Define('nBJ', 'Sum(maskBJClean)')\
            .Filter('nJClean >= 4 && nBJ > 0 && nBJ <= 2', 'GoodJet')
            
    return fdf
#.Filter('nCleanJet >= 4 && nBJet > 0 && nBJet <= 2', 'GoodJet')\

# - iSkim3 ########################################################################
def FSkim3(df):
    fdf = df.Filter('iSkim == 3', 'iSkim3')\
            .Define('maskMu', 'Muon_pfRelIso03_all < 0.15 && Muon_pt > 15')\
            .Define('nGoodMu', 'Sum(maskMu)')\
            .Define('maskEl', 'Electron_pfRelIso03_all < 0.15')\
            .Define('nGoodEl', 'Sum(maskEl)')\
            .Filter('nGoodMu == 3 && nGoodEl == 0', 'GoodEvent')\
            .Filter('abs(Sum(Muon_charge[maskMu])) != 3', 'GoodCharge')\
            .Define('Muon_pt15', 'Muon_pt[maskMu]')\
            .Define('maskJClean', 'cleanJ(Jet_eta, Muon_eta[maskMu], Jet_phi, Muon_phi[maskMu])')\
            .Define('maskBjet', 'Jet_btagDeepB > 0.6')\
            .Define('nJClean', 'Sum(maskJClean)')\
            .Define('maskBJClean', 'prod(maskJClean,maskBjet)')\
            .Define('nBJ', 'Sum(maskBJClean)')\
            .Filter('nJClean >= 2 && nBJ >= 1', 'GoodJet')
    
    return fdf
#.Filter('nCleanJet >= 2 && nBJet >= 1', 'GoodJet')\

# - iSkim4 ########################################################################
def FSkim4(df):
    fdf = df.Filter('iSkim == 4', 'iSkim4')\
            .Define('maskMu', 'Muon_pfRelIso03_all < 0.15 && Muon_pt > 15')\
            .Define('nGoodMu', 'Sum(maskMu)')\
            .Define('Muon_pt15', 'Muon_pt[maskMu]')\
            .Define('maskEl', 'Electron_pfRelIso03_all < 0.15 && Electron_pt > 15')\
            .Define('nGoodEl', 'Sum(maskEl)')\
            .Define('Electron_pt15', 'Electron_pt[maskEl]')\
            .Filter('nGoodMu == 2 && nGoodEl == 1', 'GoodEvent')\
            .Define('ChargeMus', 'Muon_charge[maskMu][0] * Muon_charge[maskMu][1]')\
            .Filter('ChargeMus < 0', 'GoodCharge')\
            .Define('maskJClean', 'cleanJ(Jet_eta, Merge(Muon_eta[maskMu], Electron_eta[maskEl]),\
            Jet_phi,Merge(Muon_phi[maskMu], Electron_phi[maskEl]))')\
            .Define('maskBjet', 'Jet_btagDeepB > 0.7')\
            .Define('nJClean', 'Sum(maskJClean)')\
            .Define('maskBJClean', 'prod(maskJClean,maskBjet)')\
            .Define('nBJ', 'Sum(maskBJClean)')\
            .Filter('nJClean >= 2 && nBJ >= 1', 'GoodJet')
            
    return fdf
#.Filter('nCleanJet >= 2 && nBJet >= 1', 'GoodJet')\

# - iSkim dictionary ##############################################################

FSkim ={1 : FSkim1, 2 : FSkim2, 3 : FSkim3, 4 : FSkim4}

# - Declare Variables ####################################################################

def DeclareVariables1(df, title, save=True):
    finalVariables1 = {'mu_pt0','mu_pt1','mu_mass0','mu_mass1','max_eta','mu_q0','mu_q1','inv_m',
                       'MET_pt','nJClean','nBJ','dRHJ','dRLJ','ST','dPhi0','dPhi1','eventWeightLumi'}
    
    define =  FSkim1(df).Define('mu_pt0', 'Muon_pt[maskMu][0]')\
                        .Define('mu_pt1', 'Muon_pt[maskMu][1]')\
                        .Define('mu_eta0', 'Muon_eta[maskMu][0]')\
                        .Define('mu_eta1', 'Muon_eta[maskMu][1]')\
                        .Define('max_eta', 'max(abs(mu_eta0), abs(mu_eta1))')\
                        .Define('mu_phi0', 'Muon_phi[maskMu][0]')\
                        .Define('mu_phi1', 'Muon_phi[maskMu][1]')\
                        .Define('mu_mass0', '0.1057')\
                        .Define('mu_mass1', '0.1057')\
                        .Define('mu_q0', 'Muon_charge[maskMu][0]')\
                        .Define('mu_q1', 'Muon_charge[maskMu][1]')\
                        .Define('inv_m', 'InvMass2(mu_pt0, mu_pt1, mu_eta0, mu_eta1, mu_phi0, mu_phi1, mu_mass0, mu_mass1)')\
                        .Filter('inv_m > 15', 'GoodMass')\
                        .Filter('abs(inv_m - 91.2) > 10 && mu_q0 * mu_q1 < 0', 'rmZ')\
                        .Define('dRHJ', 'closest_Jet(Jet_eta[maskJClean],mu_eta0,Jet_phi[maskJClean],mu_phi0)')\
                        .Define('dRLJ', 'closest_Jet(Jet_eta[maskJClean],mu_eta1,Jet_phi[maskJClean],mu_phi1)')\
                        .Define('ST', 'mu_pt0 + mu_pt1 + Sum(Jet_pt[maskJClean]) + MET_pt')\
                        .Define('dPhi0', 'dPhi(MET_phi, mu_phi0)')\
                        .Define('dPhi1', 'dPhi(MET_phi, mu_phi1)')\
                        .Define('notBJet','prod(maskJClean, (-1*maskBJClean)+1)')\
                        .Define('di_Jet_invm', 'diJ_invm(Jet_pt[notBJet], Jet_eta[notBJet], Jet_phi[notBJet],\
                                Jet_mass[notBJet])')
    
    if save: define.Snapshot('Events', dirOutPath + title + 'Flat1.root', finalVariables1)
    
    return define

##############################################################################################

def DeclareVariables2(df, title, save=True):
    finalVariables2 = {'mu_pt','el_pt','mu_mass','el_mass','max_eta','mu_q','el_q','inv_m',
                       'MET_pt','nJClean', 'nBJ', 'dRHJ', 'dRLJ', 'ST','dPhi0','dPhi1','eventWeightLumi'}

    define =  FSkim2(df).Define('mu_pt', 'Muon_pt[maskMu][0]')\
                        .Define('el_pt', 'Electron_pt[maskEl][0]')\
                        .Define('mu_eta', 'Muon_eta[maskMu][0]')\
                        .Define('el_eta', 'Electron_eta[maskEl][0]')\
                        .Define('max_eta', 'max(abs(mu_eta), abs(el_eta))')\
                        .Define('mu_phi', 'Muon_phi[maskMu][0]')\
                        .Define('el_phi', 'Electron_phi[maskEl][0]')\
                        .Define('mu_mass', '0.1057')\
                        .Define('el_mass', '0.000511')\
                        .Define('mu_q', 'Muon_charge[maskMu][0]')\
                        .Define('el_q', 'Electron_charge[maskEl][0]')\
                        .Define('inv_m', 'InvMass2(mu_pt, el_pt, mu_eta, el_eta, mu_phi, el_phi, mu_mass, el_mass)')\
                        .Filter('inv_m > 15', 'GoodMass')\
                        .Define('maxELepton', 'higherELepton(Merge(Muon_pt[maskMu],Electron_pt[maskEl]),Merge(Muon_eta[maskMu],\
                                Electron_eta[maskEl]),Merge(Muon_phi[maskMu], Electron_phi[maskEl]))')\
                        .Define('minELepton', 'higherELepton(Merge(Muon_pt[maskMu],Electron_pt[maskEl]),Merge(Muon_eta[maskMu],\
                                Electron_eta[maskEl]),Merge(Muon_phi[maskMu], Electron_phi[maskEl]), 1)')\
                        .Define('dRHJ', 'closest_Jet(Jet_eta[maskJClean],maxELepton[0],Jet_phi[maskJClean],maxELepton[1])')\
                        .Define('dRLJ', 'closest_Jet(Jet_eta[maskJClean],minELepton[0],Jet_phi[maskJClean],minELepton[1])')\
                        .Define('ST', 'mu_pt + el_pt + Sum(Jet_pt[maskJClean]) + MET_pt')\
                        .Define('dPhi0', 'dPhi(MET_phi, mu_phi)')\
                        .Define('dPhi1', 'dPhi(MET_phi, el_phi)')\
                        .Define('notBJet','prod(maskJClean, (-1*maskBJClean)+1)')\
                        .Define('di_Jet_invm', 'diJ_invm(Jet_pt[notBJet], Jet_eta[notBJet], Jet_phi[notBJet],\
                                Jet_mass[notBJet])')
                        
    if save: define.Snapshot('Events', dirOutPath + title + 'Flat2.root', finalVariables2)
    
    return define

##############################################################################################

def DeclareVariables3(df, title, save=True):
    finalVariables3 = {'mu_pt0','mu_pt1','mu_pt2','mu_mass0','mu_mass1','mu_mass2','mu_q0','mu_q1','mu_q2','inv_m01',
                       'inv_m12','inv_m02','inv_m3','MET_pt','nJClean','nBJ','Jet_pt0','Jet_pt1','bJet_pt0','dR1J',
                       'dR0bJ','dR01','dR02','ST','dPhi0','dPhi1','dPhi2','eventWeightLumi'}

    define =  FSkim3(df).Define('mu_idx', 'find_idx(Muon_charge[maskMu], Muon_phi[maskMu], Muon_eta[maskMu])')\
                        .Define('mu_pt0', 'Muon_pt[maskMu][mu_idx[0]]')\
                        .Define('mu_pt1', 'Muon_pt[maskMu][mu_idx[1]]')\
                        .Define('mu_pt2', 'Muon_pt[maskMu][mu_idx[2]]')\
                        .Define('mu_eta0', 'Muon_eta[maskMu][mu_idx[0]]')\
                        .Define('mu_eta1', 'Muon_eta[maskMu][mu_idx[1]]')\
                        .Define('mu_eta2', 'Muon_eta[maskMu][mu_idx[2]]')\
                        .Define('mu_phi0', 'Muon_phi[maskMu][mu_idx[0]]')\
                        .Define('mu_phi1', 'Muon_phi[maskMu][mu_idx[1]]')\
                        .Define('mu_phi2', 'Muon_phi[maskMu][mu_idx[2]]')\
                        .Define('mu_mass0', '0.1057')\
                        .Define('mu_mass1', '0.1057')\
                        .Define('mu_mass2', '0.1057')\
                        .Define('mu_q0', 'Muon_charge[maskMu][mu_idx[0]]')\
                        .Define('mu_q1', 'Muon_charge[maskMu][mu_idx[1]]')\
                        .Define('mu_q2', 'Muon_charge[maskMu][mu_idx[2]]')\
                        .Define('inv_m01', 'InvMass2(mu_pt0, mu_pt1, mu_eta0, mu_eta1, mu_phi0, mu_phi1, mu_mass0, mu_mass1)')\
                        .Define('inv_m12', 'InvMass2(mu_pt1, mu_pt2, mu_eta1, mu_eta2, mu_phi1, mu_phi2, mu_mass1, mu_mass2)')\
                        .Define('inv_m02', 'InvMass2(mu_pt0, mu_pt2, mu_eta0, mu_eta2, mu_phi0, mu_phi2, mu_mass0, mu_mass2)')\
                        .Define('inv_m3', 'InvMass3(Muon_pt[maskMu], Muon_eta[maskMu], Muon_phi[maskMu], Muon_mass[maskMu])')\
                        .Filter('inv_m01 > 15')\
                        .Filter('inv_m02 > 15')\
                        .Filter('inv_m12 > 15')\
                        .Filter('abs(inv_m01 - 91.2) > 10 && abs(inv_m02 - 91.2) > 10', 'rmZ')\
                        .Define('dR01', 'dR(mu_eta0, mu_eta1, mu_phi0, mu_phi1)')\
                        .Define('dR02', 'dR(mu_eta0, mu_eta2, mu_phi0, mu_phi2)')\
                        .Define('Jet_pt0', 'Jet_pt[maskJClean][0]')\
                        .Define('Jet_pt1', 'Jet_pt[maskJClean][1]')\
                        .Define('bJet_pt0', 'Jet_pt[maskBJClean][0]')\
                        .Define('dR1J', 'closest_Jet(Jet_eta[maskJClean],mu_eta1,Jet_phi[maskJClean],mu_phi1)')\
                        .Define('dR0bJ', 'closest_Jet(Jet_eta[maskBJClean],mu_eta0,Jet_phi[maskBJClean],mu_phi0)')\
                        .Define('ST', 'mu_pt0 + mu_pt1  + mu_pt2 + Sum(Jet_pt[maskJClean]) + MET_pt')\
                        .Define('dPhi0', 'dPhi(MET_phi, mu_phi0)')\
                        .Define('dPhi1', 'dPhi(MET_phi, mu_phi1)')\
                        .Define('dPhi2', 'dPhi(MET_phi, mu_phi2)')\
                        .Define('notBJet','prod(maskJClean, (-1*maskBJClean)+1)')\
                        .Define('di_Jet_invm', 'diJ_invm(Jet_pt[notBJet], Jet_eta[notBJet], Jet_phi[notBJet],\
                                Jet_mass[notBJet])')
    
    if save: define.Snapshot('Events', dirOutPath + title + 'Flat3.root', finalVariables3)
    
    return define

#.Filter('abs(inv_m02 - 91.2) > 10', 'other Z contamination with one lepton not reconstructed')\
#.Define('WJ','abs(di_Jet_invm - 80) < 5')\
#.Filter('Sum(WJ) < 1', 'cutting W')

##############################################################################################

def DeclareVariables4(df, title, save=True):    
    finalVariables4 = {'lep_pt0','lep_pt1','lep_pt2','lep_mass0','lep_mass1','lep_mass2','lep_q0','lep_q1','lep_q2','inv_m01',
                       'inv_m12','inv_m02','inv_m3','MET_pt','nJClean','nBJ','Jet_pt0','Jet_pt1','bJet_pt0','dR1J',
                       'dR0bJ','dR01','dR02','ST','dPhi0','dPhi1','dPhi2','eventWeightLumi'}

    define =  FSkim4(df).Define('lep_pt', 'Merge(Muon_pt[maskMu], Electron_pt[maskEl])')\
                        .Define('lep_eta', 'Merge(Muon_eta[maskMu], Electron_eta[maskEl])')\
                        .Define('lep_phi', 'Merge(Muon_phi[maskMu], Electron_phi[maskEl])')\
                        .Define('lep_mass', 'RVec<float> {0.1057, 0.1057, 0.000511}')\
                        .Define('lep_charge', 'Merge(Muon_charge[maskMu], Electron_charge[maskEl])')\
                        .Define('lep_idx', 'find_idx(lep_charge, lep_phi, lep_eta)')\
                        .Define('lep_pt0', 'lep_pt[lep_idx[0]]')\
                        .Define('lep_pt1', 'lep_pt[lep_idx[1]]')\
                        .Define('lep_pt2', 'lep_pt[lep_idx[2]]')\
                        .Define('lep_eta0', 'lep_eta[lep_idx[0]]')\
                        .Define('lep_eta1', 'lep_eta[lep_idx[1]]')\
                        .Define('lep_eta2', 'lep_eta[lep_idx[2]]')\
                        .Define('lep_phi0', 'lep_phi[lep_idx[0]]')\
                        .Define('lep_phi1', 'lep_phi[lep_idx[1]]')\
                        .Define('lep_phi2', 'lep_phi[lep_idx[2]]')\
                        .Define('lep_mass0', 'lep_mass[lep_idx[0]]')\
                        .Define('lep_mass1', 'lep_mass[lep_idx[1]]')\
                        .Define('lep_mass2', 'lep_mass[lep_idx[2]]')\
                        .Define('lep_q0', 'lep_charge[lep_idx[0]]')\
                        .Define('lep_q1', 'lep_charge[lep_idx[1]]')\
                        .Define('lep_q2', 'lep_charge[lep_idx[2]]')\
                        .Define('inv_m01', 'InvMass2(lep_pt0, lep_pt1, lep_eta0, lep_eta1, lep_phi0, lep_phi1, lep_mass0,\
                        lep_mass1)')\
                        .Define('inv_m12', 'InvMass2(lep_pt1, lep_pt2, lep_eta1, lep_eta2, lep_phi1, lep_phi2, lep_mass1,\
                        lep_mass2)')\
                        .Define('inv_m02', 'InvMass2(lep_pt0, lep_pt2, lep_eta0, lep_eta2, lep_phi0, lep_phi2, lep_mass0,\
                        lep_mass2)')\
                        .Define('inv_m3', 'InvMass3(lep_pt, lep_eta, lep_phi, lep_mass)')\
                        .Filter('abs(inv_m01 - 91.2) > 10 && abs(inv_m02 - 91.2) > 10', 'rmZ')\
                        .Define('dR01', 'dR(lep_eta0, lep_eta1, lep_phi0, lep_phi1)')\
                        .Define('dR02', 'dR(lep_eta0, lep_eta2, lep_phi0, lep_phi2)')\
                        .Define('Jet_pt0', 'Jet_pt[maskJClean][0]')\
                        .Define('Jet_pt1', 'Jet_pt[maskJClean][1]')\
                        .Define('bJet_pt0', 'Jet_pt[maskBJClean][0]')\
                        .Define('dR1J', 'closest_Jet(Jet_eta[maskJClean],lep_eta1,Jet_phi[maskJClean],lep_phi1)')\
                        .Define('dR0bJ', 'closest_Jet(Jet_eta[maskBJClean],lep_eta0,Jet_phi[maskBJClean],lep_phi0)')\
                        .Define('ST', 'lep_pt0 + lep_pt1  + lep_pt2 + Sum(Jet_pt[maskJClean]) + MET_pt')\
                        .Define('dPhi0', 'dPhi(MET_phi, lep_phi0)')\
                        .Define('dPhi1', 'dPhi(MET_phi, lep_phi1)')\
                        .Define('dPhi2', 'dPhi(MET_phi, lep_phi2)')\
                        .Define('notBJet','prod(maskJClean, (-1*maskBJClean)+1)')\
                        .Define('di_Jet_invm', 'diJ_invm(Jet_pt[notBJet], Jet_eta[notBJet], Jet_phi[notBJet],\
                                Jet_mass[notBJet])')
    
    if save: define.Snapshot('Events', dirOutPath + title + 'Flat4.root', finalVariables4)
    
    return define

##############################################################################################

DeclareVariables = {1 : DeclareVariables1, 2 : DeclareVariables2, 3 : DeclareVariables3, 4 : DeclareVariables4}