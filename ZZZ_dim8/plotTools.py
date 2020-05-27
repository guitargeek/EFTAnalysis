from LHEevent import *
import ROOT as rt

def InvariantMassSq(p1, p2):
    v1 = rt.TLorentzVector(p1['Px'], p1['Py'], p1['Pz'], p1['E'])
    v2 = rt.TLorentzVector(p2['Px'], p2['Py'], p2['Pz'], p2['E'])
    return (v1+v2).Mag2()
