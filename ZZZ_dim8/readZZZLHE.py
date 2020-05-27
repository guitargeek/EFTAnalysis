import sys
import ROOT as rt
import math
from LHEevent import *
from LHEfile import *
import plotTools

if __name__ == '__main__':

    #Bprime histograms
    h_Lepton_1_Pt = rt.TH1F("h_Lepton_1_Pt", "h_Lepton_1_Pt", 1000, 0., 1000)
    h_Lepton_1_Pt.Sumw2()
    h_Lepton_2_Pt = rt.TH1F("h_Lepton_2_Pt", "h_Lepton_2_Pt", 1000, 0., 1000)
    h_Lepton_2_Pt.Sumw2()
    h_Lepton_3_Pt = rt.TH1F("h_Lepton_3_Pt", "h_Lepton_3_Pt", 1000, 0., 1000)
    h_Lepton_3_Pt.Sumw2()
    h_Lepton_4_Pt = rt.TH1F("h_Lepton_4_Pt", "h_Lepton_4_Pt", 1000, 0., 1000)
    h_Lepton_4_Pt.Sumw2()
    h_Lepton_5_Pt = rt.TH1F("h_Lepton_5_Pt", "h_Lepton_5_Pt", 1000, 0., 1000)
    h_Lepton_5_Pt.Sumw2()
    h_Lepton_6_Pt = rt.TH1F("h_Lepton_6_Pt", "h_Lepton_6_Pt", 1000, 0., 1000)
    h_Lepton_6_Pt.Sumw2()
    h_DRll = rt.TH1F("h_DRll", "h_DRll", 500, 0., 5.0)
    h_DRll.Sumw2()
    lep1_lv = rt.TLorentzVector()
    lep2_lv = rt.TLorentzVector()
    lep3_lv = rt.TLorentzVector()
    lep4_lv = rt.TLorentzVector()
    lep5_lv = rt.TLorentzVector()
    lep6_lv = rt.TLorentzVector()
    h_M6l = rt.TH1F("h_M6l", "h_M6l,", 5000, 0., 5000)
    h_M6l.Sumw2()
    cross = rt.TVector3()
    Z1_lv = rt.TVector3()
    Z2_lv = rt.TVector3()
    Z3_lv = rt.TVector3()
    h_acop = rt.TH1F("h_acop", "h_acop", 200, 0.0, 1.0)

    # find events in file
    myLHEfile = LHEfile(sys.argv[1])
    myLHEfile.setMax(100000)
    #myLHEfile.setMax(2)
    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        n_l = 0
        n_Z = 0
        transverseMomentum = []
        lep_lv = []
        deltaRTmp = 1000.0
        DRll = -999
        DR = []
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            if abs(p['ID']) == 23:
               n_Z += 1
               if (n_Z==1): normOne = rt.TMath.Sqrt(pow(p['Px'], 2) + pow(p['Py'], 2) + pow(p['Pz'], 2))
               if (n_Z==1): Z1_lv = rt.TVector3(p['Px']/normOne, p['Py']/normOne, p['Pz']/normOne)
               if (n_Z==2): normTwo = rt.TMath.Sqrt(pow(p['Px'], 2) + pow(p['Py'], 2) + pow(p['Pz'], 2))
               if (n_Z==2): Z2_lv = rt.TVector3(p['Px']/normTwo, p['Py']/normTwo, p['Pz']/normTwo)
               if (n_Z==3): normThree = rt.TMath.Sqrt(pow(p['Px'], 2) + pow(p['Py'], 2) + pow(p['Pz'], 2))
               if (n_Z==3): Z3_lv = rt.TVector3(p['Px']/normThree, p['Py']/normThree, p['Pz']/normThree)
               if (n_Z==3): cross = Z1_lv.Cross(Z2_lv)
               if (n_Z==3): acop =  Z3_lv.Dot(cross)
               if (n_Z==3): print acop
               if (n_Z==3): h_acop.Fill(acop)
            if (abs(p['ID']) == 11 or abs(p['ID']) == 13 or abs(p['ID']) == 15): 
                transverseMomentum.append(rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py']))
                transverseMomentum.sort(reverse=True)
                #print transverseMomentum
                pT = rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py'])
                n_l += 1
                #if n_l==6:  print n_l
                if n_l==1:  lep1_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
                if n_l==2:  lep2_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
                if n_l==3:  lep3_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
                if n_l==4:  lep4_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
                if n_l==5:  lep5_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
                if n_l==6:  lep6_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
                if n_l==6:  h_Lepton_1_Pt.Fill(transverseMomentum[0])
                if n_l==6:  h_Lepton_2_Pt.Fill(transverseMomentum[1])
                if n_l==6:  h_Lepton_3_Pt.Fill(transverseMomentum[2])
                if n_l==6:  h_Lepton_4_Pt.Fill(transverseMomentum[3])
                if n_l==6:  h_Lepton_5_Pt.Fill(transverseMomentum[4])
                if n_l==6:  h_Lepton_6_Pt.Fill(transverseMomentum[5])
                if n_l==6:  h_M6l.Fill((lep1_lv+lep2_lv+lep3_lv+lep4_lv+lep5_lv+lep6_lv).M())
                lep_lv.append(rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
                #if(len(lep_lv) > 0 and n_l==6): print lep_lv[0].Pt()
                #if(len(lep_lv) > 0 and n_l==6): h_Mass_6l.Fill()
                #print len(lep_lv)
                for i in range(0, len(lep_lv)):
                    for j in range(i+1, len(lep_lv)):
                        #print i, j
                        deltaR = lep_lv[i].DeltaR(lep_lv[j])
                        if(deltaR < deltaRTmp):
                            deltaRTmp = deltaR
                            DRll = deltaRTmp
                            DR.append(DRll)
                            DR.sort(reverse=False)
                            #print DRll
                if(n_l==6): h_DRll.Fill(DRll);
                #if(n_l==6): print DRll, DR[0]

        del oneEvent, myLHEevent
        
    # write the histograms
    histoFILE = rt.TFile(sys.argv[2],"RECREATE")
    h_Lepton_1_Pt.Write()
    h_Lepton_2_Pt.Write()
    h_Lepton_3_Pt.Write()
    h_Lepton_4_Pt.Write()
    h_Lepton_5_Pt.Write()
    h_Lepton_6_Pt.Write()
    h_DRll.Write()
    h_M6l.Write()
    h_acop.Write()
    histoFILE.Close()
