## ############################################################### ##
## RC_FIAP (Reinforced Concrete Frame Inelastic Analysis Platform) ##
##                                                                 ##
## Developed by:                                                   ##
##       Victor F. Ceballos (vceballos@uninorte.edu.co)            ##
##       Carlos A. Arteta (carteta@uninorte.edu.co                 ##
## RCFIAPMain.py : this is the main script that calls              ##
## GUIFrameNonLinearACI.py : graphical environment                 ##
## mplwidget.py : cript to help plot the plastic hinge projector   ##
## ############################################################### ##

import sys
from math import pi, sqrt, ceil, floor
from scipy import interpolate

import openseespy.opensees as op
from PyQt5.QtWidgets import *
# from PyQt5.uic import loadUi
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt

from PyQt5.QtWidgets import QDialog, QApplication
from GUIFrameNonLinearACI1d import *
import numpy as np  # load the numpy module, calling it np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
matplotlib.use("Qt5Agg")
import pandas as pd
import os
import subprocess
import runpy

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Definition of units
m = 1.  # define basic units -- output units
kN = 1.  # define basic units -- output units
sec = 1.  # define basic units -- output units
mm = m / 1000.  # define engineering units
cm = m / 100.
N = kN / 1000.
MPa = N / mm ** 2
GPa = MPa * 1000
m2 = m ** 2  # m^2
m3 = m ** 3  # m^3
m4 = m ** 4  # m^4
inch = cm * 2.54
ft = 12. * inch
g = 9.81 * m / sec ** 2  # gravitational acceleration
kip = 4.448 * kN
ksi = kip / inch ** 2
psi = ksi / 1000.
lbf = psi * inch ** 2  # pounds force
pcf = lbf / ft ** 3  # pounds per cubic foot
psf = lbf / ft ** 3  # pounds per square foot
in2 = inch ** 2  # inch^2
in4 = inch ** 4  # inch^4
GConc = 24. * kN / m ** 3  # Specific gravity of concrete
cbar = False
np.set_printoptions(precision=6)

class RegistroBeams:
    def __init__(self, tbl_data_design_beams, id_, b, h, L_As_top, L_As_bot, L_Leg_n, L_Sstirrup, R_As_top, R_As_bot,
                 R_Leg_n, R_Sstirrup):
        fila = tbl_data_design_beams.rowCount()
        tbl_data_design_beams.insertRow(fila)

        self.spx_id = QLineEdit(tbl_data_design_beams) # setattr(self, 'spx_id', QLineEdit(tbl_data_design_beams))
        self.spx_id.setValidator(QIntValidator(0, 100))
        self.spx_id.setText(f'B{id_}')
        # self.spx_id.setStyleSheet('border-top: none; border-right: none; border-bottom: none')
        # self.spx_id.setFont(('Times', 10))

        self.spx_b = QLineEdit(tbl_data_design_beams)
        self.spx_b.setValidator(QIntValidator(20, 1000))
        self.spx_b.setText('{:d}'.format(int(b)))
        # self.spx_b.setStyleSheet('border-top: none; border-right: none; border-bottom: none')

        self.spx_h = QLineEdit(tbl_data_design_beams)
        self.spx_h.setValidator(QIntValidator(20, 1000))
        self.spx_h.setText('{:d}'.format(int(h)))
        # self.spx_h.setStyleSheet('border-top: none; border-right: none; border-bottom: none')

        self.spx_L_As_top = QLineEdit(tbl_data_design_beams)
        self.spx_L_As_top.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_L_As_top.setText('{:.2f}'.format(L_As_top))
        # self.spx_L_As_top.setStyleSheet('border-top: none; border-right: none; border-bottom: none')

        self.spx_L_As_bot = QLineEdit(tbl_data_design_beams)
        self.spx_L_As_bot.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_L_As_bot.setText('{:.2f}'.format(L_As_bot))

        self.spx_L_Leg_n = QLineEdit(tbl_data_design_beams)
        self.spx_L_Leg_n.setValidator(QIntValidator(2, 10))
        self.spx_L_Leg_n.setText('{:d}'.format(int(L_Leg_n)))

        self.spx_L_Sstirrup = QLineEdit(tbl_data_design_beams)
        self.spx_L_Sstirrup.setValidator(QIntValidator(4, 30))
        self.spx_L_Sstirrup.setText('{:d}'.format(int(L_Sstirrup)))

        self.spx_R_As_top = QLineEdit(tbl_data_design_beams)
        self.spx_R_As_top.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_R_As_top.setText('{:.2f}'.format(R_As_top))

        self.spx_R_As_bot = QLineEdit(tbl_data_design_beams)
        self.spx_R_As_bot.setValidator(QDoubleValidator(2., 400., 2))
        self.spx_R_As_bot.setText('{:.2f}'.format(R_As_bot))

        self.spx_R_Leg_n = QLineEdit(tbl_data_design_beams)
        self.spx_R_Leg_n.setValidator(QIntValidator(2, 10))
        self.spx_R_Leg_n.setText('{:d}'.format(int(R_Leg_n)))

        self.spx_R_Sstirrup = QLineEdit(tbl_data_design_beams)
        self.spx_R_Sstirrup.setValidator(QIntValidator(4, 30))
        self.spx_R_Sstirrup.setText('{:d}'.format(int(R_Sstirrup)))

        tbl_data_design_beams.setCellWidget(fila, 0, self.spx_id)
        tbl_data_design_beams.setCellWidget(fila, 1, self.spx_b)
        tbl_data_design_beams.setCellWidget(fila, 2, self.spx_h)

        tbl_data_design_beams.setCellWidget(fila, 3, self.spx_L_As_top)
        tbl_data_design_beams.setCellWidget(fila, 4, self.spx_L_As_bot)
        tbl_data_design_beams.setCellWidget(fila, 5, self.spx_L_Leg_n)
        tbl_data_design_beams.setCellWidget(fila, 6, self.spx_L_Sstirrup)

        tbl_data_design_beams.setCellWidget(fila, 7, self.spx_R_As_top)
        tbl_data_design_beams.setCellWidget(fila, 8, self.spx_R_As_bot)
        tbl_data_design_beams.setCellWidget(fila, 9, self.spx_R_Leg_n)
        tbl_data_design_beams.setCellWidget(fila, 10, self.spx_R_Sstirrup)

        tbl_data_design_beams.setColumnWidth(0, 40)
        tbl_data_design_beams.setColumnWidth(1, 40)
        tbl_data_design_beams.setColumnWidth(2, 40)
        tbl_data_design_beams.setColumnWidth(3, 60)
        tbl_data_design_beams.setColumnWidth(4, 60)
        tbl_data_design_beams.setColumnWidth(5, 60)
        tbl_data_design_beams.setColumnWidth(6, 60)
        tbl_data_design_beams.setColumnWidth(7, 60)
        tbl_data_design_beams.setColumnWidth(8, 60)
        tbl_data_design_beams.setColumnWidth(9, 60)
        tbl_data_design_beams.setColumnWidth(10, 60)


class RegistroColumns:
    def __init__(self, tbl_data_design_columns, id_, b, h, db, nbH, nbB, Leg_n_H, Leg_n_B, Sstirrup):
        fila = tbl_data_design_columns.rowCount()
        tbl_data_design_columns.insertRow(fila)

        self.spx_id = QLineEdit(tbl_data_design_columns)
        self.spx_id.setValidator(QIntValidator(0, 1000))
        self.spx_id.setText(f'C{id_}')

        self.spx_b = QLineEdit(tbl_data_design_columns)
        self.spx_b.setValidator(QIntValidator(20, 1000))
        self.spx_b.setText('{:d}'.format(int(b)))

        self.spx_h = QLineEdit(tbl_data_design_columns)
        self.spx_h.setValidator(QIntValidator(20, 1000))
        self.spx_h.setText('{:d}'.format(int(h)))

        self.spx_db = QLineEdit(tbl_data_design_columns)
        self.spx_db.setValidator(QDoubleValidator(1., 10., 2))
        self.spx_db.setText('{:.2f}'.format(db))

        self.spx_nbH = QLineEdit(tbl_data_design_columns)
        self.spx_nbH.setValidator(QIntValidator(2, 100))
        self.spx_nbH.setText('{:d}'.format(int(nbH)))

        self.spx_nbB = QLineEdit(tbl_data_design_columns)
        self.spx_nbB.setValidator(QIntValidator(2, 100))
        self.spx_nbB.setText('{:d}'.format(int(nbB)))

        self.spx_Leg_n_H = QLineEdit(tbl_data_design_columns)
        self.spx_Leg_n_H.setValidator(QIntValidator(2, 100))
        self.spx_Leg_n_H.setText('{:d}'.format(int(Leg_n_H)))

        self.spx_Leg_n_B = QLineEdit(tbl_data_design_columns)
        self.spx_Leg_n_B.setValidator(QIntValidator(2, 100))
        self.spx_Leg_n_B.setText('{:d}'.format(int(Leg_n_B)))

        self.spx_Sstirrup = QLineEdit(tbl_data_design_columns)
        self.spx_Sstirrup.setValidator(QIntValidator(2, 100))
        self.spx_Sstirrup.setText('{:d}'.format(int(Sstirrup)))

        tbl_data_design_columns.setCellWidget(fila, 0, self.spx_id)
        tbl_data_design_columns.setCellWidget(fila, 1, self.spx_b)
        tbl_data_design_columns.setCellWidget(fila, 2, self.spx_h)
        tbl_data_design_columns.setCellWidget(fila, 3, self.spx_db)
        tbl_data_design_columns.setCellWidget(fila, 4, self.spx_nbH)
        tbl_data_design_columns.setCellWidget(fila, 5, self.spx_nbB)
        tbl_data_design_columns.setCellWidget(fila, 6, self.spx_Leg_n_H)
        tbl_data_design_columns.setCellWidget(fila, 7, self.spx_Leg_n_B)
        tbl_data_design_columns.setCellWidget(fila, 8, self.spx_Sstirrup)

        tbl_data_design_columns.setColumnWidth(0, 40)
        tbl_data_design_columns.setColumnWidth(1, 40)
        tbl_data_design_columns.setColumnWidth(2, 40)
        tbl_data_design_columns.setColumnWidth(3, 60)
        tbl_data_design_columns.setColumnWidth(4, 40)
        tbl_data_design_columns.setColumnWidth(5, 40)
        tbl_data_design_columns.setColumnWidth(6, 60)
        tbl_data_design_columns.setColumnWidth(7, 60)
        tbl_data_design_columns.setColumnWidth(8, 60)


class BeamElasticElement:
    def __init__(self, EleTag, Nod_ini, Nod_end, AEle, EcEle, IzEle, LEle, BEle, HEle, ElegTr, RZi, RZe):
        self.EleTag = EleTag
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.AEle = AEle
        self.EcEle = EcEle
        self.IzEle = IzEle
        self.LEle = LEle
        self.BEle = BEle
        self.HEle = HEle
        self.ElegTr = ElegTr
        self.RZi = RZi
        self.RZe = RZe


class BeamDesing:
    def __init__(self, EleTag, b, h, Ast1, dt1, Mn_n1, Asb1, db1, Mn_p1, ns1, ss1, Ast2, dt2, Mn_n2, Asb2, db2, Mn_p2,
                 ns2, ss2, Nod_ini, Nod_end, db_t1, db_b1, db_t2, db_b2, Vpr, VU1, VU2):
        self.EleTag = EleTag
        self.b = b
        self.h = h
        self.Ast1 = Ast1
        self.dt1 = dt1
        self.Mn_n1 = Mn_n1
        self.Asb1 = Asb1
        self.db1 = db1
        self.Mn_p1 = Mn_p1
        self.ns1 = ns1
        self.ss1 = ss1
        self.Ast2 = Ast2
        self.dt2 = dt2
        self.Mn_n2 = Mn_n2
        self.Asb2 = Asb2
        self.db2 = db2
        self.Mn_p2 = Mn_p2
        self.ns2 = ns2
        self.ss2 = ss2
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.db_t1 = db_t1
        self.db_b1 = db_b1
        self.db_t2 = db_t2
        self.db_b2 = db_b2
        self.Vpr = Vpr
        self.VU1 = VU1
        self.VU2 = VU2


class ColDesing:
    def __init__(self, EleTag, b, h, nbH, nbB, db, As, Pu_v, Mu_v, fiPn, fiMn, Mn_i, d, dist, ro, Mu_i,
                 sst, nsB, nsH, Nod_ini, Nod_end, NUD1, NUD2, NUG1, NUG2, MUD1, MUD2, VUD1, VUD2, ColBeamStr):
        self.EleTag = EleTag
        self.b = b
        self.h = h
        self.nbH = nbH
        self.nbB = nbB
        self.db = db
        self.As = As
        self.Pu_v = Pu_v
        self.Mu_v = Mu_v
        self.fiPn = fiPn
        self.fiMn = fiMn
        self.Mn_i = Mn_i
        self.d = d
        self.dist = dist
        self.ro = ro
        self.Mu_i = Mu_i
        self.sst = sst
        self.nsB = nsB
        self.nsH = nsH
        self.Nod_ini = Nod_ini
        self.Nod_end = Nod_end
        self.NUD1 = NUD1
        self.NUD2 = NUD2
        self.NUG1 = NUG1
        self.NUG2 = NUG2
        self.MUD1 = MUD1
        self.MUD2 = MUD2
        self.VUD1 = VUD1
        self.VUD2 = VUD2
        self.ColBeamStr = ColBeamStr


class DuctilityCurve:
    def __init__(self, xi, xe, yi, ye, CD_i, CD_e):
        self.xi = xi
        self.xe = xe
        self.yi = yi
        self.ye = ye
        self.CD_i = CD_i
        self.CD_e = CD_e


class PlasticRotationAngle:
    def __init__(self, xi, xe, yi, ye, PRA_i, PRA_e):
        self.xi = xi
        self.xe = xe
        self.yi = yi
        self.ye = ye
        self.PRA_i = PRA_i
        self.PRA_e = PRA_e


class AcceptanceCriteria:
    def __init__(self, IO_1, LS_1, CP_1, IO_2, LS_2, CP_2):
        self.IO_1 = IO_1
        self.LS_1 = LS_1
        self.CP_1 = CP_1
        self.IO_2 = IO_2
        self.LS_2 = LS_2
        self.CP_2 = CP_2


class PlasticHingeLength:
    def __init__(self, phl1, phl2):
        self.phl1 = phl1
        self.phl2 = phl2


class MyForm(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_NonLinearFrameDialog()
        #self.setStyleSheet("QLineEdit {border: none}")
        self.ui.setupUi(self)
        self.ui.Design.clicked.connect(self.Design)
        self.ui.CreateNLM.clicked.connect(self.CreateNLM)
        self.ui.Pushover.clicked.connect(self.Pushover)
        self.ui.IDA.clicked.connect(self.IDA)
        self.ui.progressBarPushover.hide()
        self.ui.progressBarIDA.hide()
        self.ui.progressBarBeamDesign.hide()
        self.ui.progressBarColumnDesign.hide()
        self.ui.Exit.clicked.connect(self.Exit)
        self.show()

    def Exit(self):
        self.close()

    def Design(self):
        # global Loc_span, Loc_heigth, ListNodes, Elements, DataBeamDesing, DataColDesing, WDL, WLL, WDLS, Wtotal, \
        #     cover, num_elems, Beta1B, Beta1C, fcB, fcC
        CodeDesign = self.ui.comboBoxDesignCode.currentText()
        if CodeDesign == 'ACI 318S-19 IMF':
            exec(open("Design_ACI_318S_19_IFM.py").read())
        elif CodeDesign == 'ACI 318S-19 SMF':
            exec(open("Design_ACI_318S_19_SFM.py").read())

    # Creation of the nonlinear model
    def CreateNLM(self):
        global T1m, T2m, EleCol, EleBeam, MG_ElemsForceS1, MG_ElemsDeforS1, MG_ElemsForceS6, MG_ElemsDeforS6,\
            DataBeamPhl, DataColPhl

        # Validation of beam and column design table data
        def validate_data(self):
            cover = 4*cm
            dst = 3/8*inch
            for (r, DB) in zip(self.registros_beams, DataBeamDesing):
                DB.b = float(r.spx_b.text()) * cm
                DB.h = float(r.spx_h.text()) * cm
                DB.Ast1 = float(r.spx_L_As_top.text()) * cm ** 2
                DB.Asb1 = float(r.spx_L_As_bot.text()) * cm ** 2
                DB.Ast2 = float(r.spx_R_As_top.text()) * cm ** 2
                DB.Asb2 = float(r.spx_R_As_bot.text()) * cm ** 2
                DB.ns1 = int(r.spx_L_Leg_n.text())
                DB.ns2 = int(r.spx_R_Leg_n.text())
                DB.ss1 = float(r.spx_L_Sstirrup.text()) * cm
                DB.ss2 = float(r.spx_R_Sstirrup.text()) * cm

            for (r, DC) in zip(self.registros_cols, DataColDesing):
                DC.b = float(r.spx_b.text()) * cm
                DC.h = float(r.spx_h.text()) * cm
                DC.db = float(r.spx_db.text()) * mm
                DC.nbH = int(r.spx_nbH.text())
                DC.nbB = int(r.spx_nbB.text())
                DC.nsH = int(r.spx_Leg_n_H.text())
                DC.nsB = int(r.spx_Leg_n_B.text())
                DC.sst = float(r.spx_Sstirrup.text()) * cm
                dp = cover + dst + 0.5 * DC.db
                DC.dist = np.linspace(dp, DC.h - dp, DC.nbH)
                Ab = pi * DC.db ** 2. / 4.
                DC.As = np.hstack([DC.nbB * Ab, np.ones(DC.nbH - 2) * 2 * Ab, DC.nbB * Ab])

        # Function: Parameters of regularized unconfined concrete
        def con_inconf_regu():
            fpc = -fc
            epsc0 = 2 * fpc / Ec
            Gfc = max(2.0 * (-fpc / MPa) * N / mm, 25.0 * N / mm)
            epscu = Gfc / (0.6 * fpc * phl) - 0.8 * fpc / Ec + epsc0
            fcu = 0.2 * fpc
            lambdaU = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            return fpc, epsc0, fcu, epscu, lambdaU, ft, Ets

        # Function: Parameters of regularized confined concrete
        def con_conf_regu(b, h, nsB, nsH, sst):
            fpc = -fc
            bcx = h - 2. * cover - dst
            bcy = b - 2. * cover - dst
            Asx = nsB * Ast
            Asy = nsH * Ast
            Asvt = Asx + Asy
            flx = Asvt * fy / sst / bcx
            fly = Asvt * fy / sst / bcy
            slx = bcx / (nsB - 1)
            sly = bcy / (nsH - 1)
            k2x = min(0.26 * sqrt((bcx / sst) * (bcx / slx) * (1000. / flx)), 1)
            k2y = min(0.26 * sqrt((bcy / sst) * (bcy / sly) * (1000. / fly)), 1)
            flex = k2x * flx
            fley = k2y * fly
            fle = (flex * bcx + fley * bcy) / (bcx + bcy)
            k1 = 6.7 * (fle / 1000.) ** (-0.17)
            fcc = fc + k1 * fle
            fpcc = -fcc
            Ecc = Ec
            Gfc = max(2.0 * (-fpc / MPa) * N / mm, 25.0 * N / mm)
            K = k1 * fle / fc
            epscc0 = eo1 * (1. + 5. * K)
            Gfcc = 1.7 * Gfc
            epsccu = Gfcc / (0.6 * fpcc * phl) - 0.8 * fpcc / Ecc + epscc0
            fccu = 0.2 * fpcc
            lambdaC = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            return fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets

        # Function: Regularized steel parameters
        def steel_mat_regu():
            FyTestN4 = 490.0 * MPa
            FsuTestN4 = 630.0 * MPa
            epsuTestN4 = 0.10
            LgageTestN4 = 200.0 * mm
            Es = 200.0 * GPa
            FyPosN4 = FyTestN4
            epsyPosN4 = FyPosN4 / Es
            FyNegN4 = FyTestN4
            epsyNegN4 = FyNegN4 / Es
            FsuPosN4 = FsuTestN4
            epsuPosN4 = epsyPosN4 + LgageTestN4 / phl * (epsuTestN4 - epsyPosN4)
            bPosN4 = (FsuPosN4 - FyPosN4) / (Es * (epsuPosN4 - epsyPosN4))
            epsuNegN4 = min(-epsccu, epsuPosN4)
            bNegN4 = bPosN4
            # FsuNegN4 = FsuTestN4
            FsuNegN4 = FyNegN4 + bNegN4 * (Es * (epsuNegN4 - epsyNegN4))
            FsrPosN4 = 0.2 * FyPosN4
            epsrPosN4 = (FsuPosN4 - FsrPosN4) / Es + 1.05 * epsuPosN4
            FsrNegN4 = 0.2 * FsuNegN4
            epsrNegN4 = (FsuNegN4 - FsrNegN4) / Es + 1.05 * epsuNegN4
            pinchX = 0.2
            pinchY = 0.8
            damage1 = 0.0
            damage2 = 0.0
            beta = 0.0
            # op.uniaxialMaterial('Hysteretic', Ele.EleTag * 6 + 4 + pos, FyPosN4, epsyPosN4, FsuPosN4, epsuPosN4
            #                     , FsrPosN4, epsrPosN4, -FyNegN4, -epsyNegN4, -FsuNegN4, -epsuNegN4, -FsrNegN4
            #                     , -epsrNegN4, pinchX, pinchY, damage1, damage2, beta)

            SteelN4Mat = Ele.EleTag * 6 + 4 + pos
            SteelMPFTag = int(1e6 * SteelN4Mat)
            R0 = 20.0
            cR1 = 0.925
            cR2 = 0.15
            a1 = 0.0
            a2 = 1.0
            a3 = 0.0
            a4 = 0.0
            print(Ele.EleTag, 'SteelMPF', int(SteelMPFTag), FyPosN4/MPa, FyNegN4/MPa, Es/GPa, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3, a4)
            op.uniaxialMaterial('SteelMPF', SteelMPFTag, FyPosN4, FyNegN4, Es, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3,
                                a4)
            print('MinMax', int(SteelN4Mat), int(SteelMPFTag), '-min', -epsuNegN4, '-max', epsuPosN4)
            op.uniaxialMaterial('MinMax', SteelN4Mat, SteelMPFTag, '-min', -epsuNegN4, '-max', epsuPosN4)

        # Function: Parameters of non-regularized confined concrete
        def con_conf(b, h, nsB, nsH, sst):
            fpc = -fc
            bcx = h - 2. * cover - dst
            bcy = b - 2. * cover - dst
            Asx = nsB * Ast
            Asy = nsH * Ast
            Asvt = Asx + Asy
            flx = Asvt * fy / sst / bcx
            fly = Asvt * fy / sst / bcy
            slx = bcx / (nsB - 1)
            sly = bcy / (nsH - 1)
            k2x = min(0.26 * sqrt((bcx / sst) * (bcx / slx) * (1000. / flx)), 1)
            k2y = min(0.26 * sqrt((bcy / sst) * (bcy / sly) * (1000. / fly)), 1)
            flex = k2x * flx
            fley = k2y * fly
            fle = (flex * bcx + fley * bcy) / (bcx + bcy)
            k1 = 6.7 * (fle / 1000.) ** (-0.17)
            fcc = fc + k1 * fle
            fpcc = -fcc
            K = k1 * fle / fc
            epscc0 = eo1 * (1. + 5. * K)
            rov = Asvt / sst / (bcx + bcy)
            e85 = 260 * rov * epscc0 + eo85
            epsccu = (e85 - epscc0) * (0.2 * fcc - fcc) / (0.85 * fcc - fcc) + epscc0
            fccu = 0.2 * fpcc
            lambdaC = 0.10
            ft = 0.33 * sqrt(-fpc * MPa)
            Ets = ft / 0.002
            return fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets

        # Function: Parameters of non-regularized steel
        def steel_mat():
            FyTestN4 = 490.0 * MPa
            FsuTestN4 = 630.0 * MPa
            epsuTestN4 = 0.10
            LgageTestN4 = phl
            Es = 200.0 * GPa
            FyPosN4 = FyTestN4
            epsyPosN4 = FyPosN4 / Es
            FyNegN4 = FyTestN4
            epsyNegN4 = FyNegN4 / Es
            FsuPosN4 = FsuTestN4
            epsuPosN4 = epsyPosN4 + LgageTestN4 / phl * (epsuTestN4 - epsyPosN4)
            bPosN4 = (FsuPosN4 - FyPosN4) / (Es * (epsuPosN4 - epsyPosN4))
            epsuNegN4 = min(-epsccu, epsuPosN4)
            bNegN4 = bPosN4
            # FsuNegN4 = FsuTestN4
            FsuNegN4 = FyNegN4 + bNegN4 * (Es * (epsuNegN4 - epsyNegN4))
            FsrPosN4 = 0.2 * FyPosN4
            epsrPosN4 = (FsuPosN4 - FsrPosN4) / Es + 1.05 * epsuPosN4
            FsrNegN4 = 0.2 * FsuNegN4
            epsrNegN4 = (FsuNegN4 - FsrNegN4) / Es + 1.05 * epsuNegN4
            pinchX = 0.2
            pinchY = 0.8
            damage1 = 0.0
            damage2 = 0.0
            beta = 0.0
            # op.uniaxialMaterial('Hysteretic', Ele.EleTag * 6 + 4 + pos, FyPosN4, epsyPosN4, FsuPosN4, epsuPosN4
            #                     , FsrPosN4, epsrPosN4, -FyNegN4, -epsyNegN4, -FsuNegN4, -epsuNegN4, -FsrNegN4
            #                     , -epsrNegN4, pinchX, pinchY, damage1, damage2, beta)

            SteelN4Mat = Ele.EleTag * 6 + 4 + pos
            SteelMPFTag = 1e6 * SteelN4Mat
            R0 = 20.0
            cR1 = 0.925
            cR2 = 0.15
            a1 = 0.0
            a2 = 1.0
            a3 = 0.0
            a4 = 0.0
            # print('SteelMPF', int(SteelMPFTag), FyPosN4, FyNegN4, Es, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3, a4)
            op.uniaxialMaterial('SteelMPF', SteelMPFTag, FyPosN4, FyNegN4, Es, bPosN4, bNegN4, R0, cR1, cR2, a1, a2, a3,
                                a4)
            # print('MinMax', int(SteelN4Mat), int(SteelMPFTag), '-min', -epsuNegN4, '-max', epsuPosN4)
            op.uniaxialMaterial('MinMax', SteelN4Mat, SteelMPFTag, '-min', -epsuNegN4, '-max', epsuPosN4)

        # Function: Creation of fibers in beams
        def fiber_beam(Ast, Asb, pos):
            op.section('Fiber', Ele.EleTag * 2 + pos)
            op.patch('rect', Ele.EleTag * 6 + 2 + pos, 10, 1, -y2 + dp, -z2 + dp, y2 - dp, z2 - dp)
            op.patch('rect', Ele.EleTag * 6 + pos, 10, 1, -y2 + dp, z2 - dp, y2 - dp, z2)
            op.patch('rect', Ele.EleTag * 6 + pos, 10, 1, -y2 + dp, -z2, y2 - dp, -z2 + dp)
            op.patch('rect', Ele.EleTag * 6 + pos, 2, 1, -y2, -z2, -y2 + dp, z2)
            op.patch('rect', Ele.EleTag * 6 + pos, 2, 1, y2 - dp, -z2, y2, z2)
            print(Ele.EleTag * 6 + 4 + pos, 1, Ast/cm**2, y2 - dp, z2 - dp, y2 - dp, -z2 + dp)
            op.layer('straight', Ele.EleTag * 6 + 4 + pos, 1, Ast, y2 - dp, z2 - dp, y2 - dp, -z2 + dp)
            print(Ele.EleTag * 6 + 4 + pos, 1, Asb/cm**2, -y2 + dp, z2 - dp, -y2 + dp, -z2 + dp)
            op.layer('straight', Ele.EleTag * 6 + 4 + pos, 1, Asb, -y2 + dp, z2 - dp, -y2 + dp, -z2 + dp)


        validate_data(self)
        op.wipe()  # The models is restarted in opensees
        op.model('Basic', '-ndm', 2, '-ndf', 3)
        for node in ListNodes:
            op.node(int(node[0]), int(node[1]), int(node[2]))
            if node[2] == 0.:
                op.fix(int(node[0]), 1, 1, 1)
        cover = 4 * cm
        dst = 3 / 8 * inch
        Ast = pi * dst ** 2 / 4.  # area de la barra del estribo
        # creacion de columnas
        HBeam = float(self.ui.HBeam.text())
        HColi = float(self.ui.HColi.text())  # Column inside Depth            # if node[2] > 0 and node[1] == 0:
            #     MasterNode = node[0]
            # if node[2] > 0 and node[1] != 0:
            #     op.equalDOF(int(MasterNode), int(node[0]), 1)

        HCole = float(self.ui.HCole.text())  # Column outside Depth
        fy = float(self.ui.fy.text()) * MPa
        Es = 200.0 * GPa
        fcB = float(self.ui.fcB.text()) * MPa
        fcC = float(self.ui.fcC.text()) * MPa
        op.geomTransf('PDelta', 1, '-jntOffset', 0, 0, 0, -HBeam / 2)
        op.geomTransf('PDelta', 2, '-jntOffset', 0, HBeam / 2, 0, -HBeam / 2)
        op.geomTransf('Corotational', 3, '-jntOffset', HColi / 2., 0, -HColi / 2., 0)
        op.geomTransf('Corotational', 4, '-jntOffset', HCole / 2., 0, -HColi / 2., 0)
        op.geomTransf('Corotational', 5, '-jntOffset', HColi / 2., 0, -HCole / 2., 0)

        EleCol = []
        EleBeam = []
        for Ele in Elements:
            if ListNodes[Ele.Nod_ini, 1] == ListNodes[Ele.Nod_end, 1]:
                EleCol.append(Ele)
            else:
                EleBeam.append(Ele)

        # Creation of non-linear elements (beams and columns)
        eo1, eo85, eo20, lambdaU = -0.002, -0.0038, -0.006, 0.1
        DataColPhl = []
        for (Ele, DC) in zip(EleCol, DataColDesing):
            fc, Ec = fcC, Ele.EcEle
            if self.ui.radioButton05H.isChecked() == True:
                phl = 0.5 * DC.h
            if self.ui.radioButtonPark.isChecked() == True:
                phl = 0.08 * Ele.LEle + 0.022 * fy / MPa * DC.db / mm
            if self.ui.radioButtonBerry.isChecked() == True:
                phl = 0.05 * Ele.LEle + 0.1 * fy / MPa * DC.db / mm / sqrt(fc * MPa)
            DataColPhl.append(PlasticHingeLength(phl, phl))
            if self.ui.radioButtonYesRegu.isChecked() == True:
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                print('fpc, epsc0, fcu, epscu, lambdaU, ft, Ets', fpc/MPa, epsc0, fcu/MPa, epscu, lambdaU, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 1, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_regu(DC.b, DC.h, DC.nsB, DC.nsH, DC.sst)
                print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc/MPa, epscc0, fccu/MPa, epsccu, lambdaC, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                pos = 0
                steel_mat_regu()
                pos = 1
                steel_mat_regu()
            if self.ui.radioButtonNoRegu.isChecked() == True:
                ft = 0.33 * sqrt(fcC * MPa)
                Ets = ft / 0.002
                # print('Concrete02', Ele.EleTag * 6, -fcC, eo1, -0.2 * fcC, eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6, -fcC, eo1, -0.2 * fcC, eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 1, -fcC, eo1, -0.2 * fcC, eo20, lambdaU, ft, Ets)
                fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DC.b, DC.h, DC.nsB, DC.nsH, DC.sst)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                pos = 0
                steel_mat()
                pos = 1
                steel_mat()

            dp = DC.dist[0]
            y1 = DC.h / 2.0
            z1 = DC.b / 2.0
            op.section('Fiber', Ele.EleTag)
            op.patch('rect', Ele.EleTag * 6 + 2, 10, 1, -y1 + dp, -z1 + dp, y1 - dp, z1 - dp)
            op.patch('rect', Ele.EleTag * 6, 10, 1, -y1 + dp, z1 - dp, y1 - dp, z1)
            op.patch('rect', Ele.EleTag * 6, 10, 1, -y1 + dp, -z1, y1 - dp, -z1 + dp)
            op.patch('rect', Ele.EleTag * 6, 2, 1, -y1, -z1, -y1 + dp, z1)
            op.patch('rect', Ele.EleTag * 6, 2, 1, y1 - dp, -z1, y1, z1)
            for dist, As in zip(DC.dist, DC.As):
                # print('Col ', Ele.EleTag * 6 + 4, 1, As, -y1 + dist, z1 - dp, -y1 + dist, -z1 + dp)
                op.layer('straight', Ele.EleTag * 6 + 4, 1, As, -y1 + dist, z1 - dp, -y1 + dist, -z1 + dp)
            MassDens = Ele.AEle * GConc / g
            op.beamIntegration('HingeRadau', Ele.EleTag, Ele.EleTag, phl, Ele.EleTag, phl, Ele.EleTag)
            op.element('forceBeamColumn', Ele.EleTag, Ele.Nod_ini, Ele.Nod_end, Ele.ElegTr, Ele.EleTag
                       , '-mass', MassDens)
        # print('DataColPhl =', DataColPhl)
        DataBeamPhl = []
        for (Ele, DB) in zip(EleBeam, DataBeamDesing):
            fc, Ec, nsH = fcB, Ele.EcEle, 2
            if self.ui.radioButton05H.isChecked() == True:
                phl1 = 0.5 * DB.h
                phl2 = 0.5 * DB.h
            if self.ui.radioButtonPark.isChecked() == True:
                phl1 = 0.08 * Ele.LEle + 0.022 * fy / MPa * DB.db_t1 / mm
                phl2 = 0.08 * Ele.LEle + 0.022 * fy / MPa * DB.db_t2 / mm
            if self.ui.radioButtonPark.isChecked() == True:
                phl1 = 0.05 * Ele.LEle + 0.1 * fy / MPa * DB.db_t1 / mm / sqrt(fc * MPa)
                phl2 = 0.05 * Ele.LEle + 0.1 * fy / MPa * DB.db_t2 / mm / sqrt(fc * MPa)
            DataBeamPhl.append(PlasticHingeLength(phl1, phl2))
            if self.ui.radioButtonYesRegu.isChecked() == True:
                phl = phl1
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                print('fpc, epsc0, fcu, epscu, lambdaU, ft, Ets', fpc/MPa, epsc0, fcu/MPa, epscu, lambdaU, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                phl = phl2
                fpc, epsc0, fcu, epscu, lambdaU, ft, Ets = con_inconf_regu()
                print('fpc, epsc0, fcu, epscu, lambdaU, ft, Ets', fpc/MPa, epsc0, fcu/MPa, epscu, lambdaU, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 1, fpc, epsc0, fcu, epscu, lambdaU, ft, Ets)
                phl, pos = phl1, 0
                fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_regu(DB.b, DB.h, DB.ns1, nsH, DB.ss1)
                print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc/MPa, epscc0, fccu/MPa, epsccu, lambdaC, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                steel_mat_regu()
                phl, pos = phl2, 1
                fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf_regu(DB.b, DB.h, DB.ns2, nsH, DB.ss2)
                print('fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets', fpcc/MPa, epscc0, fccu/MPa, epsccu, lambdaC, ft/MPa, Ets/MPa)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                steel_mat_regu()
            if self.ui.radioButtonNoRegu.isChecked() == True:
                ft = 0.33 * sqrt(fcB * MPa)
                Ets = ft / 0.002
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6, -fcB, eo1, -0.2 * fcB, eo20, lambdaU, ft, Ets)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 1, -fcB, eo1, -0.2 * fcB, eo20, lambdaU, ft, Ets)
                fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DB.b, DB.h, DB.ns1, nsH, DB.ss1)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 2, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets = con_conf(DB.b, DB.h, DB.ns2, nsH, DB.ss2)
                op.uniaxialMaterial('Concrete02', Ele.EleTag * 6 + 3, fpcc, epscc0, fccu, epsccu, lambdaC, ft, Ets)
                pos = 0
                steel_mat()
                pos = 1
                steel_mat()
            y2 = DB.h / 2.0
            z2 = DB.b / 2.0
            dp = DB.h - min(DB.db1, DB.dt1)
            print('dp1 =', dp)
            pos = 0
            fiber_beam(DB.Ast1, DB.Asb1, pos)
            dp = DB.h - min(DB.db2, DB.dt2)
            print('dp2 =', dp)
            pos = 1
            fiber_beam(DB.Ast2, DB.Asb2, pos)
            MassDens = Ele.AEle * GConc / g + WDLS / g
            print('phl1', phl1, 'phl2', phl2)
            op.beamIntegration('HingeRadau', Ele.EleTag, Ele.EleTag * 2, phl1, Ele.EleTag * 2 + 1, phl2, Ele.EleTag * 2)
            op.element('forceBeamColumn', Ele.EleTag, Ele.Nod_ini, Ele.Nod_end, Ele.ElegTr, Ele.EleTag
                       , '-mass', MassDens)

        # print('DataBeamPhl =', DataBeamPhl)
        list_beams = [Ele.EleTag for Ele in EleBeam]
        list_cols = [Ele.EleTag for Ele in EleCol]
        print('list_beams =', list_beams)
        print('list_cols =', list_cols)

        if not os.path.exists("Pushover"):
            os.mkdir("Pushover")

        # Recording of forces and deformations from nonlinear analysis
        op.recorder('Element', '-file', 'Pushover/beams_force_1.out',
                    '-time', '-ele', *list_beams, 'section', 1, 'force')
        op.recorder('Element', '-file', 'Pushover/beams_def_1.out',
                    '-time', '-ele', *list_beams, 'section', 1, 'deformation')
        op.recorder('Element', '-file', 'Pushover/beams_force_6.out',
                    '-time', '-ele', *list_beams, 'section', 6, 'force')
        op.recorder('Element', '-file', 'Pushover/beams_def_6.out',
                    '-time', '-ele', *list_beams, 'section', 6, 'deformation')
        op.recorder('Element', '-file', 'Pushover/cols_force_1.out',
                    '-time', '-ele', *list_cols, 'section', 1, 'force')
        op.recorder('Element', '-file', 'Pushover/cols_def_1.out',
                    '-time', '-ele', *list_cols, 'section', 1, 'deformation')
        op.recorder('Element', '-file', 'Pushover/cols_force_6.out',
                    '-time', '-ele', *list_cols, 'section', 6, 'force')
        op.recorder('Element', '-file', 'Pushover/cols_def_6.out',
                    '-time', '-ele', *list_cols, 'section', 6, 'deformation')
        op.recorder('Node', '-file', 'Pushover/HoriNodes.out',
                    '-time', '-node', *ListNodes, '-dof', 1, 'disp')
        op.recorder('Node', '-file', 'Pushover/VertNodes.out',
                    '-time', '-node', *ListNodes, '-dof', 2, 'disp')

        # Create a Plain load pattern for gravity loading with a Linear TimeSeries
        op.timeSeries('Linear', 1)
        op.pattern('Plain', 1, 1)
        for Ele in EleCol:
            op.eleLoad('-ele', Ele.EleTag, '-type', '-beamUniform', 0, -Ele.AEle * GConc)
        for Ele in EleBeam:
            op.eleLoad('-ele', Ele.EleTag, '-type', '-beamUniform', -Ele.AEle * GConc - WDL - 0.25*WLL)

        Tol = 1.0e-6  # convergence tolerance for test
        op.constraints('Plain')  # how it handles boundary conditions
        op.numberer('Plain')  # renumber dof to minimize band-width (optimization), if you want to
        op.system('BandGeneral')  # how to store and solve the system of equations in the analysis
        op.test('NormDispIncr', Tol, 100)  # determine if convergence has been achieved at the end of an iteration step
        op.algorithm('KrylovNewton')  # use Newton solution algorithm: updates tangent stiffness at every iteration
        NstepGravity = 10  # apply gravity in 10 steps
        DGravity = 1. / NstepGravity  # first load increment;
        op.integrator('LoadControl', DGravity)  # determine the next time step for an analysis
        op.analysis('Static')  # define type of analysis static or transient
        nele = num_elems-1
        MG_ElemsForceS1, MG_ElemsDeforS1, MG_ElemsForceS6, MG_ElemsDeforS6 = np.zeros(2*nele), np.zeros(2*nele), \
                                                                         np.zeros(2*nele), np.zeros(2*nele)
        step = 1
        loadf = 1.0

        while step <= NstepGravity and loadf > 0:
            ElemsForceS1, ElemsDeforS1, ElemsForceS6, ElemsDeforS6 = [], [], [], []
            op.analyze(1)
            for Element in Elements:
                ForcesS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'force'))
                ForcesS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'force'))
                DeforsS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'deformation'))
                DeforsS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'deformation'))
                ElemsForceS1 = np.append(ElemsForceS1, ForcesS1)
                ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
                ElemsForceS6 = np.append(ElemsForceS6, ForcesS6)
                ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
            MG_ElemsForceS1 = np.vstack((MG_ElemsForceS1, ElemsForceS1))
            MG_ElemsDeforS1 = np.vstack((MG_ElemsDeforS1, ElemsDeforS1))
            MG_ElemsForceS6 = np.vstack((MG_ElemsForceS6, ElemsForceS6))
            MG_ElemsDeforS6 = np.vstack((MG_ElemsDeforS6, ElemsDeforS6))

            loadf = op.getTime()
            step += 1
        op.loadConst('-time', 0.0)
        print("Model Nonlinear Built")

        xi = 0.02  # damping ratio
        MpropSwitch = 1.0
        KcurrSwitch = 0.0
        KcommSwitch = 1.0
        KinitSwitch = 0.0
        floors_num = len(Loc_heigth) - 1
        print('floors_num =', floors_num)
        if floors_num >= 2:
            nEigenI = 1  # mode 1
            nEigenI2 = 2  # mode 2
            nEigenJ = 3  # mode 3
            lambdaN = op.eigen(nEigenJ)  # eigenvalue analysis for nEigenJ modes
            lambdaI = lambdaN[nEigenI - 1]  # eigenvalue mode i
            lambdaI2 = lambdaN[nEigenI2 - 1]  # eigenvalue mode i2
            # lambdaJ = lambdaN[nEigenJ - 1]  # eigenvalue mode j
            print('lambdaN ', lambdaN)
            omegaI = pow(lambdaI, 0.5)
            omegaI2 = pow(lambdaI2, 0.5)
            # omegaJ = pow(lambdaJ, 0.5)
            T1m = 2. * pi / omegaI
            T2m = 2. * pi / omegaI2

            print('Ta1=', T1m, 'seg', ' Ta2=', T2m, ' seg')
            # alphaM = MpropSwitch * xi * (2. * omegaI * omegaJ) / (omegaI + omegaJ)  # M-prop. damping D = alphaM*M
            # betaKcurr = KcurrSwitch * 2. * xi / (omegaI + omegaJ)  # current-K      +beatKcurr*KCurrent
            # betaKcomm = KcommSwitch * 2. * xi / (omegaI + omegaJ)  # last-committed K   +betaKcomm*KlastCommitt
            # betaKinit = KinitSwitch * 2. * xi / (omegaI + omegaJ)  # initial-K     +beatKinit*Kini
            # op.rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm)  # RAYLEIGH damping
        elif floors_num == 1:
            nEigenI = 1  # mode 1
            # nEigenI2 = 2  # mode 2
            lambdaN = op.eigen(nEigenI)  # eigenvalue analysis for nEigenJ modes
            lambdaI = lambdaN[nEigenI - 1]  # eigenvalue mode i
            # lambdaI2 = lambdaN[nEigenI2 - 1]  # eigenvalue mode i2
            print('lambdaN ', lambdaN)
            omegaI = pow(lambdaI, 0.5)
            # omegaI2 = pow(lambdaI2, 0.5)
            T1m = 2. * pi / omegaI
            # T2m = 2. * pi / omegaI2
            print('Ta1=', T1m, 'seg')

        self.ui.tabWidget.setCurrentIndex(2)

    # Pushover function
    def Pushover(self):
        global cbar
        def singlePush1(dref, mu, ctrlNode, dispDir, nSteps):
            IOflag = 2
            testType = 'RelativeNormDispIncr'
            # set testType	EnergyIncr;					# Dont use with Penalty constraints
            # set testType	RelativeNormUnbalance;		# Dont use with Penalty constraints
            # set testType	RelativeNormDispIncr;		# Dont use with Lagrange constraints
            # set testType	RelativeTotalNormDispIncr;	# Dont use with Lagrange constraints
            # set testType	RelativeEnergyIncr;			# Dont use with Penalty constraints
            tolInit = 1.0e-6  # the initial Tolerance, so it can be referred back to
            iterInit = 50  # the initial Max Number of Iterations
            algorithmType = 'KrylovNewton'  # the algorithm type

            op.test(testType, tolInit
                    , iterInit)  # determine if convergence has been achieved at the end of an iteration step
            op.algorithm(algorithmType)  # use Newton solution algorithm: updates tangent stiffness at every iteration
            disp = dref * mu
            dU = disp / (1.0 * nSteps)
            print('dref ', dref, 'mu ', mu, 'dU ', dU, 'disp ', disp)
            op.integrator('DisplacementControl', ctrlNode, dispDir, dU)  # determine the next time step for an analysis
            op.analysis('Static')  # define type of analysis static or transient

            # Print values
            if IOflag >= 1:
                print('singlePush: Push ', ctrlNode, ' to ', mu)

            #      the initial values to start the while loop
            nele = num_elems - 1
            MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6 = np.zeros(2 * nele), np.zeros(2 * nele), \
                                                                                 np.zeros(2 * nele), np.zeros(2 * nele)
            ok = 0
            step = 1
            loadf = 1.0
            # This feature of disabling the possibility of having a negative loading has been included.
            # This has been adapted from a similar script by Prof. Garbaggio
            htot = op.nodeCoord(ctrlNode, 2)
            maxDriftPiso = 0.0
            VBasal_v = []
            DriftTecho_v = []
            while step <= nSteps and ok == 0 and loadf > 0:
                self.ui.progressBarPushover.setValue(100 * step / nSteps)
                ElemsForceS1, ElemsDeforS1, ElemsForceS6, ElemsDeforS6 = [], [], [], []
                ok = op.analyze(1)
                loadf = op.getTime()
                temp = op.nodeDisp(ctrlNode, dispDir)
                # Print the current displacement
                if IOflag >= 2:
                    print('Pushed ', ctrlNode, ' in ', dispDir, ' to ', temp, ' with ', loadf, 'step', step)

                # If the analysis fails, try the following changes to achieve convergence
                # Analysis will be slower in here though...
                if ok != 0:
                    print('Trying relaxed convergence..')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    ok = op.analyze(1)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying Newton with initial then current .')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    op.algorithm('Newton', '-initialThenCurrent')
                    ok = op.analyze(1)
                    op.algorithm(algorithmType)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying ModifiedNewton with initial ..')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    op.algorithm('ModifiedNewton', '-initial')
                    ok = op.analyze(1)
                    op.algorithm(algorithmType)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying KrylovNewton ..')
                    op.test(testType, tolInit * 0.01,
                            iterInit * 50)  # determine if convergence has been achieved at the end of an iteration step
                    op.algorithm('KrylovNewton')
                    ok = op.analyze(1)
                    op.algorithm(algorithmType)
                    op.test(testType, tolInit,
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                if ok != 0:
                    print('Trying FixedNumIter .. ....')
                    op.test('FixedNumIter',
                            iterInit)  # determine if convergence has been achieved at the end of an iteration step
                    ok = op.analyze(1)

                for (nod_ini, nod_end) in zip(ListNodesDrift[:-1, 0], ListNodesDrift[1:, 0]):
                    # print('nod_ini ', nod_ini, 'nod_end', nod_end)
                    nod_ini = int(nod_ini)
                    nod_end = int(nod_end)
                    pos_i = op.nodeCoord(nod_ini, 2)
                    pos_s = op.nodeCoord(nod_end, 2)
                    hpiso = pos_s - pos_i
                    desp_i = op.nodeDisp(nod_ini, 1)
                    desp_s = op.nodeDisp(nod_end, 1)
                    desp_piso = abs(desp_s - desp_i)
                    drift_piso = desp_piso / hpiso
                    if drift_piso >= maxDriftPiso:
                        maxDriftPiso = drift_piso

                VBasal = 0.
                op.reactions()
                for node in ListNodesBasal:
                    # print('ind Basal ', node[0])
                    VBasal = VBasal + op.nodeReaction(node[0], 1)
                VBasal_v = np.append(VBasal_v, VBasal)
                DriftTecho = op.nodeDisp(ctrlNode, dispDir) / htot
                DriftTecho_v = np.append(DriftTecho_v, DriftTecho)
                for Element in Elements:
                    ForcesS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'force'))
                    ForcesS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'force'))
                    DeforsS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'deformation'))
                    DeforsS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'deformation'))
                    ElemsForceS1 = np.append(ElemsForceS1, ForcesS1)
                    ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
                    ElemsForceS6 = np.append(ElemsForceS6, ForcesS6)
                    ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
                MP_ElemsForceS1 = np.vstack((MP_ElemsForceS1, ElemsForceS1))
                MP_ElemsDeforS1 = np.vstack((MP_ElemsDeforS1, ElemsDeforS1))
                MP_ElemsForceS6 = np.vstack((MP_ElemsForceS6, ElemsForceS6))
                MP_ElemsDeforS6 = np.vstack((MP_ElemsDeforS6, ElemsDeforS6))

                loadf = op.getTime()
                step += 1
            # print('MP_ElemsForceS1 =', MP_ElemsForceS1)
            # print('Tamñano de MP_ElemsForceS1 =', np.shape(MP_ElemsForceS1))
            maxDriftTecho = dU * step / htot
            maxDriftTecho2 = op.nodeDisp(ctrlNode, dispDir) / htot

            if ok != 0:
                print('DispControl Analysis FAILED')
            else:
                print('DispControl Analysis SUCCESSFUL')
            if loadf <= 0:
                print('Stopped because of Load factor below zero: ', loadf)
            #    if PrintFlag == 0:
            #        os.remove("singlePush.txt")
            #        print singlePush.txt
            return maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, DriftTecho_v, MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6

        # Pushover function varying tests and algorithms
        def singlePush(dref, mu, ctrlNode, dispDir, nSteps):
            # --------------------------------------------------
            # Description of Parameters
            # --------------------------------------------------
            # dref:			Reference displacement to which cycles are run. Corresponds to yield or equivalent other, such as 1mm
            # mu:			Multiple of dref to which the push is run. So pushover can be run to a specifived ductility or displacement
            # ctrlNode:		Node to control with the displacement integrator.
            # dispDir:		DOF the loading is applied.
            # nSteps:		Number of steps.
            # IOflag:		Option to print details on screen. 2 for print of each step, 1 for basic info (default), 0 for off
            # ---------------------------------------------------
            test = {1: 'NormDispIncr', 2: 'RelativeEnergyIncr', 3: 'EnergyIncr',
                    4: 'RelativeNormUnbalance', 5: 'RelativeNormDispIncr',
                    6: 'NormUnbalance', 7: 'FixedNumIter'}
            alg = {1: 'KrylovNewton', 2: 'SecantNewton', 3: 'ModifiedNewton',
                   4: 'RaphsonNewton', 5: 'PeriodicNewton', 6: 'BFGS',
                   7: 'Broyden', 8: 'NewtonLineSearch'}

            # test = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 3:'EnergyIncr'}
            # alg = {1:'KrylovNewton', 2:'ModifiedNewton'}

            IOflag = 2
            PrintFlag = 0
            testType = 'RelativeNormDispIncr'  # Dont use with Penalty constraints

            tolInit = 1.0e-7  # the initial Tolerance, so it can be referred back to
            iterInit = 50  # the initial Max Number of Iterations
            algorithmType = 'KrylovNewton'  # the algorithm type
            #      	algorithmType Newton;		#      the algorithm type
            #      	algorithmType Newton;		#      the algorithm type

            # op.constraints('Transformation') # how it handles boundary conditions
            # op.numberer('RCM')    # renumber dof to minimize band-width (optimization), if you want to
            # op.system('BandGeneral') # how to store and solve the system of equations in the analysis

            op.test(testType, tolInit,
                    iterInit)  # determine if convergence has been achieved at the end of an iteration step
            op.algorithm(algorithmType)  # use Newton solution algorithm: updates tangent stiffness at every iteration
            disp = dref * mu
            dU = disp / (1.0 * nSteps)
            print('dref ', dref, 'mu ', mu, 'dU ', dU, 'disp ', disp, 'nSteps ', nSteps)
            op.integrator('DisplacementControl', ctrlNode, dispDir, dU)  # determine the next time step for an analysis
            op.analysis('Static')  # defivne type of analysis static or transient

            # Print values
            if IOflag >= 1:
                print('singlePush: Push ', ctrlNode, ' to ', mu)
            nele = num_elems - 1
            MP_ElemsForceS1, MP_ElemsDeforS1, MP_ElemsForceS6, MP_ElemsDeforS6 = np.zeros(2 * nele), np.zeros(2 * nele), \
                                                                                 np.zeros(2 * nele), np.zeros(2 * nele)
            #      the initial values to start the while loop
            ok = 0
            step = 1
            loadf = 1.0
            # This feature of disabling the possibility of having a negative loading has been included.
            # This has been adapted from a similar script by Prof. Garbaggio
            maxDriftPiso = 0.0
            htot = op.nodeCoord(ctrlNode, 2)
            VBasal_v = []
            DriftTecho_v = []
            # factor_v = np.array([1,0.75,0.5,0.25,0.1,2,3,5,10])
            # fact_v = np.array([50,100,500])
            # factor = 100
            # fact = 1.
            while step <= nSteps and ok == 0 and loadf > 0:
                self.ui.progressBarPushover.setValue(100 * step / nSteps)
                ElemsForceS1, ElemsDeforS1, ElemsForceS6, ElemsDeforS6 = [], [], [], []
                ok = op.analyze(1)
                loadf = op.getTime()
                temp = op.nodeDisp(ctrlNode, dispDir)
                if IOflag >= 2:
                    print('Pushed ', ctrlNode, ' in ', dispDir, ' to ', temp, ' with ', loadf, 'step ', step)
                # for factor in factor_v:
                # op.integrator('DisplacementControl',ctrlNode,dispDir,factor*dU)  # determine the next time step for an analysis
                # for fact in fact_v:
                for j in alg:
                    for i in test:
                        for fact in [1, 20, 50]:
                            if ok != 0 and j >= 4 and i != 7:
                                # print('Trying ',str(alg[j]))
                                op.test(test[i], tolInit * .01, iterInit * fact)
                                op.algorithm(alg[j])
                                ok = op.analyze(1)
                                op.algorithm(algorithmType)
                                op.test(testType, tolInit, iterInit)
                            elif ok != 0 and j < 4 and i != 7:
                                # print('Trying ',str(alg[j]))
                                op.test(test[i], tolInit, iterInit * fact)
                                op.algorithm(alg[j], '-initial')
                                ok = op.analyze(1)
                                op.algorithm(algorithmType)
                                op.test(testType, tolInit, iterInit)
                            if ok == 0:
                                break
                        if ok != 0 and i == 7:
                            op.test(test[i], iterInit)
                            op.algorithm(alg[j])
                            ok = op.analyze(1)
                        if ok == 0:
                            break
                    if ok == 0:
                        break
                    # if ok == 0:
                    #     break
                    # if ok == 0:
                    #     break
                # op.integrator('DisplacementControl',ctrlNode,dispDir,dU)  # determine the next time step for an analysis
                # Calculation of maximum Drift between floors
                for (nod_ini, nod_end) in zip(ListNodesDrift[:-1, 0], ListNodesDrift[1:, 0]):
                    # print('nod_ini ', nod_ini, 'nod_end', nod_end)
                    nod_ini = int(nod_ini)
                    nod_end = int(nod_end)
                    pos_i = op.nodeCoord(nod_ini, 2)
                    pos_s = op.nodeCoord(nod_end, 2)
                    hpiso = pos_s - pos_i
                    desp_i = op.nodeDisp(nod_ini, 1)
                    desp_s = op.nodeDisp(nod_end, 1)
                    desp_piso = abs(desp_s - desp_i)
                    drift_piso = desp_piso / hpiso
                    if drift_piso >= maxDriftPiso:
                        maxDriftPiso = drift_piso

                VBasal = 0.
                op.reactions()
                for node in ListNodesBasal:
                    # print('ind Basal ', node[0])
                    VBasal = VBasal + op.nodeReaction(node[0], 1)
                VBasal_v = np.append(VBasal_v, VBasal)
                DriftTecho = op.nodeDisp(ctrlNode, dispDir) / htot
                DriftTecho_v = np.append(DriftTecho_v, DriftTecho)
                for Element in Elements:
                    ForcesS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'force'))
                    ForcesS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'force'))
                    DeforsS1 = np.array(op.eleResponse(Element.EleTag, 'section', 1, 'deformation'))
                    DeforsS6 = np.array(op.eleResponse(Element.EleTag, 'section', 6, 'deformation'))
                    ElemsForceS1 = np.append(ElemsForceS1, ForcesS1)
                    ElemsDeforS1 = np.append(ElemsDeforS1, DeforsS1)
                    ElemsForceS6 = np.append(ElemsForceS6, ForcesS6)
                    ElemsDeforS6 = np.append(ElemsDeforS6, DeforsS6)
                MP_ElemsForceS1 = np.vstack((MP_ElemsForceS1, ElemsForceS1))
                MP_ElemsDeforS1 = np.vstack((MP_ElemsDeforS1, ElemsDeforS1))
                MP_ElemsForceS6 = np.vstack((MP_ElemsForceS6, ElemsForceS6))
                MP_ElemsDeforS6 = np.vstack((MP_ElemsDeforS6, ElemsDeforS6))
                loadf = op.getTime()
                step += 1
            maxDriftTecho = dU * step / htot
            maxDriftTecho2 = op.nodeDisp(ctrlNode, dispDir) / htot

            if ok != 0:
                print('DispControl Analysis FAILED')
            else:
                print('DispControl Analysis SUCCESSFUL')
            if loadf <= 0:
                print('Stopped because of Load factor below zero: ', loadf)
            #    if PrintFlag == 0:
            #        os.remove("singlePush.txt")
            #        print singlePush.txt
            return maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, DriftTecho_v

        ListNodesDrift = ListNodes[np.where(ListNodes[:, 1] == 0.)]
        ListNodesBasal = ListNodes[np.where(ListNodes[:, 2] == 0.)]
        if T1m <= 0.5:
            k = 1.
        elif T1m <= 2.5:
            k = 0.75 + 0.5 * T1m
        else:
            k = 2.

        sumH = np.sum(np.power(Loc_heigth, k))
        floors_num = len(Loc_heigth)

        # Defining the pushover lateral distribution type
        if self.ui.radioButtonTriangular.isChecked() == True:
            Fp = np.power(Loc_heigth, k) / sumH
        if self.ui.radioButtonUniform.isChecked() == True:
            Fp = 1. / floors_num * np.ones(floors_num + 1)
        print('Fp =', Fp)
        op.loadConst('-time', 0.0)
        op.timeSeries('Linear', 2)
        op.pattern('Plain', 2, 1)
        for (node, fp, ind) in zip(ListNodesDrift, Fp, range(floors_num)):
            op.load(int(node[0]), fp, 0.0, 0.0)

        Htotal = Loc_heigth[-1]
        Der_obj = float(self.ui.Der_obj.text())
        Des_obj = Der_obj * Htotal  # Desplazamiento objetivo
        nSteps = int(self.ui.nSteps.text())
        dref = Des_obj / nSteps
        mu = nSteps
        IDctrlNode = int(ListNodesDrift[-1, 0])  # Node where displacement is read
        print('IDctrlNode =', IDctrlNode)
        IDctrlDOF = 1  # DOF x=1, y=2
        Tol = 1.0e-4  # Tolerance
        self.ui.progressBarPushover.show()
        if self.ui.radioButtonFast.isChecked() == True:
            maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, DriftTecho_v,  MP_ElemsForceS1, MP_ElemsDeforS1,\
            MP_ElemsForceS6, MP_ElemsDeforS6 = singlePush1(dref, mu, IDctrlNode, IDctrlDOF, nSteps)
        if self.ui.radioButtonForced.isChecked() == True:
            maxDriftPiso, maxDriftTecho, maxDriftTecho2, VBasal_v, DriftTecho_v,  MP_ElemsForceS1, MP_ElemsDeforS1,\
            MP_ElemsForceS6, MP_ElemsDeforS6 = singlePush(dref, mu, IDctrlNode, IDctrlDOF, nSteps)
        self.ui.progressBarPushover.hide()

        op.wipe()

        # Plot pushover curve

        fig = self.ui.PushCurve.canvas.axes
        fig.clear()
        ax = fig.add_axes([0.1, 0.2, 0.85, 0.7])
        ax.plot(DriftTecho_v*100, -VBasal_v / Wtotal, '.-')
        ax.set_ylabel('Vb/Ws')
        ax.set_xlabel('Roof Drift %')
        ax.set_title('Pushover Curve')
        ax.grid(True)
        self.ui.PushCurve.canvas.draw()
        self.ui.PushCurve.canvas.show()

        num_beams = len(EleBeam)
        num_cols = len(EleCol)
        print('num_cols =', num_cols)
        # Reading of forces and deflections of beams and columns from recorders
        M_ElemsForceS1 = np.vstack((MG_ElemsForceS1, MP_ElemsForceS1))
        M_ElemsDeforS1 = np.vstack((MG_ElemsDeforS1, MP_ElemsDeforS1))
        M_ElemsForceS6 = np.vstack((MG_ElemsForceS6, MP_ElemsForceS6))
        M_ElemsDeforS6 = np.vstack((MG_ElemsDeforS6, MP_ElemsDeforS6))

        M_BeamsForceS1 = M_ElemsForceS1[:, 2*num_cols:]
        M_BeamsDeforS1 = M_ElemsDeforS1[:, 2*num_cols:]
        M_BeamsForceS6 = M_ElemsForceS6[:, 2*num_cols:]
        M_BeamsDeforS6 = M_ElemsDeforS6[:, 2*num_cols:]
        M_ColsForceS1 = M_ElemsForceS1[:, :2*num_cols]
        M_ColsDeforS1 = M_ElemsDeforS1[:, :2*num_cols]
        M_ColsForceS6 = M_ElemsForceS6[:, :2*num_cols]
        M_ColsDeforS6 = M_ElemsDeforS6[:, :2*num_cols]
        print('M_BeamsForceS1', M_BeamsForceS1)
        print('Tamñano de M_BeamsForceS1 =', np.shape(M_BeamsForceS1))
        print('M_ColsForceS1', M_ColsForceS1)
        print('Tamñano de M_ColsForceS1 =', np.shape(M_ColsForceS1))
        # beams_force_1 = np.loadtxt('Pushover/beams_force_1.out')
        # beams_def_1 = np.loadtxt('Pushover/beams_def_1.out')
        # beams_force_6 = np.loadtxt('Pushover/beams_force_6.out')
        # beams_def_6 = np.loadtxt('Pushover/beams_def_6.out')
        # cols_force_1 = np.loadtxt('Pushover/cols_force_1.out')
        # cols_def_1 = np.loadtxt('Pushover/cols_def_1.out')
        # cols_force_6 = np.loadtxt('Pushover/cols_force_6.out')
        # cols_def_6 = np.loadtxt('Pushover/cols_def_6.out')
        # print('cols_def_1', cols_def_1)

        fy = float(self.ui.fy.text()) * MPa
        Es = 200.0 * GPa
        ey = fy / Es
        CD_Beams = np.zeros([num_beams, 2])  # Curvature Ductility - Beams
        PRA_Beams = np.zeros([num_beams, 2])  # Plastic Rotation Angle - Beams
        My_Beams = np.zeros([num_beams, 2])

        # Calculation of curvature ductility of beams and columns
        for (ind, DB, DBPhl) in zip(range(num_beams), DataBeamDesing, DataBeamPhl):
            ets_beam_1 = M_BeamsDeforS1[:-2, 2*ind] + M_BeamsDeforS1[:-2, 2*ind+1]*(DB.dt1-DB.h/2)
            ebs_beam_1 = M_BeamsDeforS1[:-2, 2*ind] + M_BeamsDeforS1[:-2, 2*ind+1]*(DB.h/2-DB.db1)
            ets_beam_6 = M_BeamsDeforS6[:-2, 2*ind] + M_BeamsDeforS6[:-2, 2*ind+1]*(DB.dt2-DB.h/2)
            ebs_beam_6 = M_BeamsDeforS6[:-2, 2*ind] + M_BeamsDeforS6[:-2, 2*ind+1]*(DB.h/2-DB.db1)
            fi_1 = np.absolute(M_BeamsDeforS1[:-2, 2 * ind + 1])
            M_beam_1 = np.absolute(M_BeamsForceS1[:-2, 2 * ind + 1])
            fi_6 = np.absolute(M_BeamsDeforS6[:-2, 2 * ind + 1])
            M_beam_6 = np.absolute(M_BeamsForceS6[:-2, 2 * ind + 1])

            # es_beam_1 = np.maximum(np.absolute(ets_beam_1), np.absolute(ebs_beam_1))
            # es_beam_6 = np.maximum(np.absolute(ets_beam_6), np.absolute(ebs_beam_6))
            # print('es_beam_1', es_beam_1, 'es_beam_6', es_beam_6)
            if np.max(ets_beam_1) <= ey and np.max(ets_beam_1) <= ey:
                CD_1 = 0
                My_1 = 0
                PRA1 = DBPhl.phl1*np.max(fi_1)
            else:
                if np.max(ets_beam_1) >= ey:
                    ft = interpolate.interp1d(ets_beam_1, M_beam_1, kind='nearest')
                    My_1t = ft(ey)
                else:
                    My_1t = float('inf')
                if np.max(ebs_beam_1) >= ey:
                    fb = interpolate.interp1d(ebs_beam_1, M_beam_1, kind='nearest')
                    My_1b = fb(ey)
                else:
                    My_1b = float('inf')
                print('ind', ind, 'My_1t', My_1t, 'My_1b', My_1b)
                My_1 = min(My_1t, My_1b)
                f = interpolate.interp1d(M_beam_1, fi_1, kind='nearest')
                fiy_1 = f(My_1)
                CD_1 = np.max(fi_1) / fiy_1
                PRA1 = DBPhl.phl1*np.max(fi_1)
            if np.max(ets_beam_6) <= ey and np.max(ebs_beam_6) <= ey:
                CD_6 = 0
                My_6 = 0
                PRA6 = DBPhl.phl2*np.max(fi_6)
            else:
                if np.max(ets_beam_6) >= ey:
                    ft = interpolate.interp1d(ets_beam_6, M_beam_6, kind='nearest')
                    My_6t = ft(ey)
                else:
                    My_6t = float('inf')
                if np.max(ebs_beam_6) >= ey:
                    fb = interpolate.interp1d(ebs_beam_6, M_beam_6, kind='nearest')
                    My_6b = fb(ey)
                else:
                    My_6b = float('inf')
                My_6 = min(My_6t, My_6b)
                f = interpolate.interp1d(M_beam_6, fi_6, kind='nearest')
                fiy_6 = f(My_6)
                CD_6 = np.max(fi_6) / fiy_6
                PRA6 = DBPhl.phl2*np.max(fi_6)
            CD_Beams[ind, :] = [CD_1, CD_6]
            PRA_Beams[ind, :] = [PRA1, PRA6]
            My_Beams[ind, :] = [My_1, My_6]

            # print('CD_Beams =', CD_Beams)

        CD_Cols = np.zeros([num_cols, 2])
        PRA_Cols = np.zeros([num_cols, 2])
        My_Cols = np.zeros([num_cols, 2])
        for (ind, DC, DCPhl) in zip(range(num_cols), DataColDesing, DataColPhl):
            ets_col_1 = np.absolute(M_ColsDeforS1[:-2, 2*ind] + M_ColsDeforS1[:-2, 2*ind+1]*(DC.d-DC.h/2))
            ebs_col_1 = np.absolute(M_ColsDeforS1[:-2, 2*ind] + M_ColsDeforS1[:-2, 2*ind+1]*(DC.h/2-DC.d))
            ets_col_6 = np.absolute(M_ColsDeforS6[:-2, 2*ind] + M_ColsDeforS6[:-2, 2*ind+1]*(DC.d-DC.h/2))
            ebs_col_6 = np.absolute(M_ColsDeforS6[:-2, 2*ind] + M_ColsDeforS6[:-2, 2*ind+1]*(DC.h/2-DC.d))
            fi_1 = np.absolute(M_ColsDeforS1[:-2, 2 * ind + 1])
            M_col_1 = np.absolute(M_ColsForceS1[:-2, 2 * ind + 1])
            fi_6 = np.absolute(M_ColsDeforS6[:-2, 2 * ind + 1])
            M_col_6 = np.absolute(M_ColsForceS6[:-2, 2 * ind + 1])
            # es_col_1 = np.maximum(np.absolute(ets_col_1), np.absolute(ebs_col_1))
            # es_col_6 = np.maximum(np.absolute(ets_col_6), np.absolute(ebs_col_6))
            # print('es_col_1', es_col_1, 'es_col_6', es_col_6)
            if np.max(ets_col_1) <= ey and np.max(ebs_col_1) <= ey:
                CD_1 = 0
                Mfy_1 = 0
                PRA1 = DCPhl.phl1*np.max(fi_1)
            else:
                if np.max(ets_col_1) >= ey:
                    ft = interpolate.interp1d(ets_col_1, M_col_1, kind='nearest')
                    Mfy_1t = ft(ey)
                else:
                    Mfy_1t = float('inf')
                if np.max(ebs_col_1) >= ey:
                    fb = interpolate.interp1d(ebs_col_1, M_col_1, kind='nearest')
                    Mfy_1b = fb(ey)
                else:
                    Mfy_1b = float('inf')
                Mfy_1 = min(Mfy_1t, Mfy_1b)
                f = interpolate.interp1d(M_col_1, fi_1, kind='nearest')
                fify_1 = f(Mfy_1)
                My_1 = np.max(M_col_1)
                fiy_1 = My_1/Mfy_1*fify_1
                CD_1 = np.max(fi_1)/fiy_1
                PRA1 = DCPhl.phl1*np.max(fi_1)
            if np.max(ets_col_6) <= ey and np.max(ebs_col_6) <= ey:
                CD_6 = 0
                Mfy_6 = 0
                PRA6 = DCPhl.phl2*np.max(fi_6)
            else:
                if np.max(ets_col_6) >= ey:
                    ft = interpolate.interp1d(ets_col_6, M_col_6, kind='nearest')
                    Mfy_6t = ft(ey)
                else:
                    Mfy_6t = float('inf')
                if np.max(ebs_col_6) >= ey:
                    fb = interpolate.interp1d(ebs_col_6, M_col_6, kind='nearest')
                    Mfy_6b = fb(ey)
                else:
                    Mfy_6b = float('inf')
                Mfy_6 = min(Mfy_6t, Mfy_6b)
                f = interpolate.interp1d(M_col_6, fi_6, kind='nearest')
                fify_6 = f(Mfy_6)
                My_6 = np.max(M_col_6)
                fiy_6 = My_6 / Mfy_6 * fify_6
                CD_6 = np.max(fi_6)/fiy_6
                PRA6 = DCPhl.phl2*np.max(fi_6)
            CD_Cols[ind, :] = [CD_1, CD_6]
            PRA_Cols[ind, :] = [PRA1, PRA6]
            My_Cols[ind, :] = [Mfy_1, Mfy_6]
            # print('CD_Cols =', CD_Cols)
        CD_Ele = np.concatenate((CD_Cols, CD_Beams), axis=0)
        PRA_Ele = np.concatenate((PRA_Cols, PRA_Beams), axis=0)
        print('CD_Ele =', CD_Ele)
        print('PRA_Ele =', PRA_Ele)

        # Drawing of curvature ductility in the plastic hinge projector
        # self.ui.PHP.canvas.figure.clear()
        fig = self.ui.PHP.canvas.axes

        if cbar:
            cbar.remove()
            fig.clear()

        # cmap = plt.cm.get_cmap("jet")
        Desp_x = np.loadtxt('Pushover/HoriNodes.out')
        Desp_y = np.loadtxt('Pushover/VertNodes.out')
        Nodes_desp_x = ListNodes[:, 1] + 3*Desp_x[-1, 1:]
        Nodes_desp_y = ListNodes[:, 2] + 3*Desp_y[-1, 1:]
        # self.figure = plt.figure()
        ax1 = fig.add_axes([0, 0.1, 1, 0.9])

        ax1.plot(Nodes_desp_x, Nodes_desp_y, 'ks')
        ax1.axis('off')

        fpos = 0.1
        fsize = 1
        DataDC = []
        DataPRA = []
        for Ele in Elements:
            xi = Nodes_desp_x[Ele.Nod_ini]
            yi = Nodes_desp_y[Ele.Nod_ini]
            xe = Nodes_desp_x[Ele.Nod_end]
            ye = Nodes_desp_y[Ele.Nod_end]
            x = np.array([xi, xe])
            y = np.array([yi, ye])
            # print(xi, yi, xe, ye)
            ax1.plot(x, y, 'k-', alpha=.3)
            Delta_x = xe-xi
            Delta_y = ye-yi
            xi_CD = xi + fpos * Delta_x
            yi_CD = yi + fpos * Delta_y
            xe_CD = xe - fpos * Delta_x
            ye_CD = ye - fpos * Delta_y
            CD_i = CD_Ele[Ele.EleTag-1, 0]
            CD_e = CD_Ele[Ele.EleTag-1, 1]
            PRA_i = PRA_Ele[Ele.EleTag-1, 0]
            PRA_e = PRA_Ele[Ele.EleTag-1, 1]
            DataDC.append(DuctilityCurve(xi_CD, xe_CD, yi_CD, ye_CD, fsize*CD_i, fsize*CD_e))
            DataPRA.append(PlasticRotationAngle(xi_CD, xe_CD, yi_CD, ye_CD, PRA_i, PRA_e))
        ax1.axis('equal')
        DC_x, DC_y, DC_size = [], [], []
        for DC in DataDC:
            DC_x.append([DC.xi, DC.xe])
            DC_y.append([DC.yi, DC.ye])
            DC_size.append([DC.CD_i, DC.CD_e])
        DC_x = np.array(DC_x)
        DC_x = DC_x.flatten()
        DC_y = np.array(DC_y)
        DC_y = DC_y.flatten()
        DC_size = np.array(DC_size)
        DC_size = DC_size.flatten()

        PRA_size = []
        for PRA in DataPRA:
            PRA_size.append([PRA.PRA_i, PRA.PRA_e])
        PRA_size = np.array(PRA_size)
        # PRA_size = DC_size.flatten()


        if self.ui.radioButtonPHP.isChecked() == True:
            cmap = plt.cm.get_cmap('jet')

            sc = ax1.scatter(DC_x, DC_y, s=3*DC_size, c=DC_size, cmap=cmap, alpha=1)

            ax2 = fig.add_axes([0.05, 0.1, 0.9, 0.03])
            norm = plt.Normalize(vmin=np.min(DC_size), vmax=np.max(DC_size))
            cbar = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal',
                                         label='Curvature ductility $\mu_{\phi}=\phi_{u}/\phi_{y}$') #, ticks=boundsTick, spacing='proportional',


            # cbar = self.ui.PHP.canvas.figure.colorbar(mappable=sc, orientation='horizontal', aspect=80) #, ticks=boundsTick, spacing='proportional',
            # cbar = self.ui.PHP.canvas.figure.colorbar(mappable=sc, orientation="horizontal", aspect=80)

            # cbar = self.ui.PHP.canvas.figure.colorbar(mappable=sc, orientation='horizontal') #, ticks=boundsTick, spacing='proportional',

            # cbar = self.ui.PHP.canvas.figure.colorbar(mappable=sc, orientation="horizontal", aspect=80)
            # ax2.set_label('Curvature ductility $\mu_{\phi}=\phi_{u}/\phi_{y}$')


            # cbar.set_label('Acceptance Criteria ASCE 41-17')
            #
            # cbar = fig.colorbar(mappable=sc, orientation="horizontal", aspect=80)
            # cbar.set_label('Curvature ductility $\mu_{\phi}=\phi_{u}/\phi_{y}$')
            # self.ui.PHP.canvas.draw()
            # self.ui.PHP.canvas.fig.add_axes([0.2, 0.2, 0.8, 0.8])  # left,bottom,width,height
            self.ui.tabWidget.setCurrentIndex(3)
            self.ui.PHP.canvas.draw()
            self.ui.PHP.canvas.show()

        if self.ui.radioButtonAC.isChecked() == True:
            # Calculation of plastic rotation angles (radians) - Performance levels
            ecu = 0.003
            ro_bal = 0.85*fcB/fy*ecu/(ecu+fy/Es)*Beta1B
            dst = 3 / 8 * inch
            Ast = pi * dst ** 2 / 4.  # area de la barra del estribo
            DataAC = []

            Cd = float(self.ui.Cd.text())
            for (ind, DC, Ele) in zip(range(num_cols), DataColDesing, EleCol):
                if Cd <= 2.0:
                    knl = 1.0
                elif Cd >= 6.0:
                    knl = 0.7
                else:
                    knl = np.interp(Cd, [2.0, 6.0], [1.0, 0.7])
                if DC.sst/DC.d <= 0.75:
                    alfa_col = 1.0
                elif DC.sst/DC.d >= 1.0:
                    alfa_col = 0.0
                else:
                    alfa_col = np.interp(DC.sst/DC.d, [0.75, 1.0], [1.0, 0.0])
                ro_t = min(max(0.075, Ast*DC.nsH/(DC.b*DC.d)), 0.0175)
                yi = ListNodes[Ele.Nod_ini, 2]
                ye = ListNodes[Ele.Nod_end, 2]
                LCol = ye - yi
                VColOE = knl*(alfa_col*Ast*DC.nsH*fy*DC.d/DC.sst +
                              0.5*sqrt(fcC*MPa)/(min(max(2, DC.MUD1/DC.VUD1), 4)/DC.d) *
                              sqrt(1 + DC.NUD1/(0.5*DC.b*DC.h*sqrt(fcC*MPa)))*0.8*DC.b*DC.h)
                Vy = 2*My_Cols[ind, 0]/LCol
                a = max(0.042 - 0.043*DC.NUD1/(DC.b*DC.h*fcC) + 0.63*ro_t - 0.023*max(Vy/VColOE, 0.2), 0)
                bc = max(0.5/(5 + DC.NUD1/(0.8*DC.b*DC.h*fcC)*1/ro_t*fcC/fy) - 0.01, a)
                if 0.1 <= DC.NUD1/(DC.b*DC.h*fcC) <= 0.5:
                    b = bc
                elif 0.5 < DC.NUD1/(DC.b*DC.h*fcC) <= 0.7:
                    b = np.interp(DC.NUD1/(DC.b*DC.h*fcC), [0.5, 0.7], [bc, max(0, a)])
                elif DC.NUD1/(DC.b*DC.h*fcC) > 0.7:
                    b = a
                elif DC.NUD1/(DC.b*DC.h*fcC) < 0.1:
                    b = max(0.5/(5 + 0.1/0.8*1/ro_t*fcC/fy) - 0.01, a)
                IO_1 = min(0.15*a, 0.005)
                LS_1 = 0.5*b
                CP_1 = 0.7*b

                VColOE = knl*(alfa_col*Ast*DC.nsH*fy*DC.d/DC.sst +
                              0.5*sqrt(fcC*MPa)/(min(max(2, DC.MUD2/DC.VUD2), 4)/DC.d) *
                              sqrt(1 + DC.NUD2/(0.5*DC.b*DC.h * sqrt(fcC*MPa)))*0.8*DC.b*DC.h)
                Vy = 2*My_Cols[ind, 1]/LCol
                a = max(0.042 - 0.043*DC.NUD2/(DC.b*DC.h*fcC) + 0.63*ro_t - 0.023*max(Vy/VColOE, 0.2), 0)
                bc = max(0.5/(5 + DC.NUD2/(0.8*DC.b*DC.h*fcC)*1/ro_t*fcC/fy) - 0.01, a)
                if 0.1 <= DC.NUD2/(DC.b*DC.h*fcC) <= 0.5:
                    b = bc
                elif 0.5 < DC.NUD2/(DC.b*DC.h*fcC) <= 0.7:
                    b = np.interp(DC.NUD1/(DC.b*DC.h*fcC), [0.5, 0.7], [bc, max(0, a)])
                elif DC.NUD2/(DC.b*DC.h*fcC) > 0.7:
                    b = a
                elif DC.NUD2/(DC.b*DC.h*fcC) < 0.1:
                    b = max(0.5/(5 + 0.1/0.8*1/ro_t*fcC/fy) - 0.01, a)
                IO_2 = min(0.15*a, 0.005)
                LS_2 = 0.5*b
                CP_2 = 0.7*b
                print('IO_1, LS_1, CP_1, IO_2, LS_2, CP_2', IO_1, LS_1, CP_1, IO_2, LS_2, CP_2)
                DataAC.append(AcceptanceCriteria(IO_1, LS_1, CP_1, IO_2, LS_2, CP_2))

            for DB in DataBeamDesing:
                ro1 = DB.Asb1/DB.b/DB.db1
                ro1p = DB.Ast1/DB.b/DB.db1
                Vs = DB.ns1*Ast*fy*DB.db1/DB.ss1
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO1 = 0.010
                    LS1 = 0.025
                    CP1 = 0.050
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO1 = 0.005
                    LS1 = 0.020
                    CP1 = 0.040
                else:
                    IO1 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    LS1 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.025, 0.020])
                    CP1 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.050, 0.040])
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO2 = 0.005
                    LS2 = 0.020
                    CP2 = 0.030
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO2 = 0.005
                    LS2 = 0.015
                    CP2 = 0.020
                else:
                    IO2 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.005])
                    LS2 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.015])
                    CP2 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.020])
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO3 = 0.005
                    LS3 = 0.020
                    CP3 = 0.030
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO3 = 0.0015
                    LS3 = 0.010
                    CP3 = 0.015
                else:
                    IO3 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS3 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.010])
                    CP3 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.015])
                if DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) <= 0.25:
                    IO4 = 0.005
                    LS4 = 0.010
                    CP4 = 0.015
                elif DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)) >= 0.50:
                    IO4 = 0.0015
                    LS4 = 0.005
                    CP4 = 0.010
                else:
                    IO4 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS4 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    CP4 = np.interp(DB.VU1 / (DB.b * DB.db1 * sqrt(fcB * MPa)), [0.25, 0.5], [0.015, 0.010])
                if DB.ss1 <= DB.db1/3 and Vs >= 3/4*DB.VU1:
                    if (ro1-ro1p)/ro_bal <= 0:
                        IO_1 = IO1
                        LS_1 = LS1
                        CP_1 = CP1
                    elif (ro1-ro1p)/ro_bal >= 0.5:
                        IO_1 = IO2
                        LS_1 = LS2
                        CP_1 = CP2
                    else:
                        IO_1 = np.interp((ro1-ro1p)/ro_bal, [0.0, 0.5], [IO1, IO2])
                        LS_1 = np.interp((ro1-ro1p)/ro_bal, [0.0, 0.5], [LS1, LS2])
                        CP_1 = np.interp((ro1-ro1p)/ro_bal, [0.0, 0.5], [CP1, CP2])
                if DB.ss1 > DB.db1/3 or Vs < 3/4*DB.VU1:
                    if (ro1-ro1p)/ro_bal <= 0:
                        IO_1 = IO3
                        LS_1 = LS3
                        CP_1 = CP3
                    elif (ro1-ro1p)/ro_bal >= 0.5:
                        IO_1 = IO4
                        LS_1 = LS4
                        CP_1 = CP4
                    else:
                        IO_1 = np.interp((ro1-ro1p)/ro_bal, [0.0, 0.5], [IO3, IO4])
                        LS_1 = np.interp((ro1-ro1p)/ro_bal, [0.0, 0.5], [LS3, LS4])
                        CP_1 = np.interp((ro1-ro1p)/ro_bal, [0.0, 0.5], [CP3, CP4])
                ro2 = DB.Ast2 / DB.b / DB.db2
                ro2p = DB.Asb2 / DB.b / DB.db2
                Vs = DB.ns2 * Ast * fy * DB.db2 / DB.ss2
                if DB.VU2 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO1 = 0.010
                    LS1 = 0.025
                    CP1 = 0.050
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO1 = 0.005
                    LS1 = 0.020
                    CP1 = 0.040
                else:
                    IO1 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    LS1 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.025, 0.020])
                    CP1 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.050, 0.040])
                if DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO2 = 0.005
                    LS2 = 0.020
                    CP2 = 0.030
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO2 = 0.005
                    LS2 = 0.015
                    CP2 = 0.020
                else:
                    IO2 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.005])
                    LS2 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.015])
                    CP2 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.020])
                if DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO3 = 0.005
                    LS3 = 0.020
                    CP3 = 0.030
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO3 = 0.0015
                    LS3 = 0.010
                    CP3 = 0.015
                else:
                    IO3 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS3 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.020, 0.010])
                    CP3 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.030, 0.015])
                if DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) <= 0.25:
                    IO4 = 0.005
                    LS4 = 0.010
                    CP4 = 0.015
                elif DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)) >= 0.50:
                    IO4 = 0.0015
                    LS4 = 0.005
                    CP4 = 0.010
                else:
                    IO4 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.005, 0.0015])
                    LS4 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.010, 0.005])
                    CP4 = np.interp(DB.VU1 / (DB.b * DB.db2 * sqrt(fcB * MPa)), [0.25, 0.5], [0.015, 0.010])
                if DB.ss2 <= DB.db2 / 3 and Vs >= 3 / 4 * DB.VU2:
                    if (ro2 - ro2p) / ro_bal <= 0:
                        IO_2 = IO1
                        LS_2 = LS1
                        CP_2 = CP1
                    elif (ro2 - ro2p) / ro_bal >= 0.5:
                        IO_2 = IO2
                        LS_2 = LS2
                        CP_2 = CP2
                    else:
                        IO_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [IO1, IO2])
                        LS_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [LS1, LS2])
                        CP_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [CP1, CP2])
                if DB.ss1 > DB.db2 / 3 or Vs < 3 / 4 * DB.VU1:
                    if (ro2 - ro2p) / ro_bal <= 0:
                        IO_2 = IO3
                        LS_2 = LS3
                        CP_2 = CP3
                    elif (ro2 - ro2p) / ro_bal >= 0.5:
                        IO_2 = IO4
                        LS_2 = LS4
                        CP_2 = CP4
                    else:
                        IO_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [IO3, IO4])
                        LS_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [LS3, LS4])
                        CP_2 = np.interp((ro2 - ro2p) / ro_bal, [0.0, 0.5], [CP3, CP4])
                print('IO_1, LS_1, CP_1, IO_2, LS_2, CP_2', IO_1, LS_1, CP_1, IO_2, LS_2, CP_2)
                DataAC.append(AcceptanceCriteria(IO_1, LS_1, CP_1, IO_2, LS_2, CP_2))

            # print('DataAC =\n', DataAC)
            AC_size = []
            for AC in DataAC:
                AC_size.append([AC.IO_1, AC.LS_1, AC.CP_1, AC.IO_2, AC.LS_2, AC.CP_2])
            AC_size = np.array(AC_size)
            print('AC_size', AC_size)
            Color_PRA = []
            for ind in range(num_elems - 1):
                print('ind =', ind)
                if PRA_size[ind, 0] <= AC_size[ind, 0]:
                    Color_PRA.append('green')
                elif AC_size[ind, 0] < PRA_size[ind, 0] <= AC_size[ind, 1]:
                    Color_PRA.append('yellow')
                elif AC_size[ind, 1] < PRA_size[ind, 0] <= AC_size[ind, 2]:
                    Color_PRA.append('orange')
                elif PRA_size[ind, 0] > AC_size[ind, 2]:
                    Color_PRA.append('red')
                if PRA_size[ind, 1] <= AC_size[ind, 3]:
                    Color_PRA.append('green')
                elif AC_size[ind, 3] < PRA_size[ind, 1] <= AC_size[ind, 4]:
                    Color_PRA.append('yellow')
                elif AC_size[ind, 4] < PRA_size[ind, 1] <= AC_size[ind, 5]:
                    Color_PRA.append('orange')
                elif PRA_size[ind, 1] > AC_size[ind, 5]:
                    Color_PRA.append('red')
            # cm = plt.cm.get_cmap('jet')

            from matplotlib.colors import ListedColormap
            # cMap = ListedColormap(['green', 'yellow', 'orange', 'red'])

            cmap = colors.ListedColormap(['green', 'yellow', 'orange', 'red'])
            norm = colors.BoundaryNorm([0, 1, 2, 3, 4], cmap.N)

            PRA_size = DC_size.flatten()
            ax1.scatter(DC_x, DC_y, s=30, c=Color_PRA, cmap=cmap, norm=norm, alpha=1)
            # divider = make_axes_locatable(self.ui.PHP.canvas.axes)
            # cax = divider.new_vertical(size="5%", pad=0.5, pack_start=True)
            ax2 = fig.add_axes([0.05, 0.1, 0.9, 0.03])
            # cmap = colors.ListedColormap(['green', 'yellow', 'orange', 'red'])
            # norm = colors.BoundaryNorm([0, 1, 2, 3, 4], cmap.N)
            cmap.set_over('0.25')
            cmap.set_under('0.75')
            bounds = [0, 1, 2, 3, 4]
            boundsTick = ['0', 'IO', 'LS', 'CP']
            norm = colors.BoundaryNorm(bounds, cmap.N)
            cbar = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation='horizontal') #, ticks=boundsTick, spacing='proportional',
                                        # )
            cbar.set_label('Acceptance Criteria ASCE 41-17')
            #
            # cbar = self.ui.PHP.canvas.figure.colorbar(mappable=sc, orientation="horizontal") #, aspect=80)
            # cbar.set_label('Curvature ductility $\mu_{\phi}=\phi_{u}/\phi_{y}$')
            # # cbar.set_ticks([0, 1, 2, 3, 4])
            cbar.set_ticklabels(['', 'IO', 'LS', 'CP', ''])
            self.ui.PHP.canvas.draw()
            self.ui.PHP.canvas.show()

            # self.ui.PHP.canvas.axes.set_position([0.2, 0.3, 0.9, 0.9])  # left,bottom,width,height
            self.ui.tabWidget.setCurrentIndex(3)

            # # The second example illustrates the use of a ListedColormap, a
            # # BoundaryNorm, and extended ends to show the "over" and "under"
            # # value colors.
            # cmap = mpl.colors.ListedColormap(['r', 'g', 'b', 'c'])
            # cmap.set_over('0.25')
            # cmap.set_under('0.75')
            #
            # # If a ListedColormap is used, the length of the bounds array must be
            # # one greater than the length of the color list.  The bounds must be
            # # monotonically increasing.
            # bounds = [1, 2, 4, 7, 8]
            # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            # cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
            #                                 norm=norm,
            #                                 # to use 'extend', you must
            #                                 # specify two extra boundaries:
            #                                 boundaries=[0] + bounds + [13],
            #                                 extend='both',
            #                                 ticks=bounds,  # optional
            #                                 spacing='proportional',
            #                                 orientation='horizontal')
            # cb2.set_label('Discrete intervals, some other units')

            # cMap = ListedColormap(['white', 'green', 'blue', 'red'])
            #
            # # data
            # np.random.seed(42)
            # data = np.random.rand(4, 4)
            # fig, ax = plt.subplots()
            # heatmap = ax.pcolor(data, cmap=cMap)
            #
            # # legend
            # cbar = plt.colorbar(heatmap)
            #
            # cbar.ax.get_yaxis().set_ticks([])
            # for j, lab in enumerate(['$0$', '$1$', '$2$', '$>3$']):
            #     cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='center', va='center')
            # cbar.ax.get_yaxis().labelpad = 15
            # cbar.ax.set_ylabel('# of contacts', rotation=270)

        # Plot of the lower moment-curvature diagram of the first column
        indb = 1 #+ 2 * (naxes - 1)
        indc = 2

        fig = self.ui.MomFi.canvas.axes
        fig.clear()
        ax = fig.add_axes([0.1, 0.2, 0.85, 0.7])
        ax.plot(np.absolute(M_ColsDeforS6[:-2, 2*(indc-1)+1]), np.absolute(M_ColsForceS6[:-2, 2*(indc-1)+1]), 'b.-')
        ax.plot(np.absolute(M_BeamsDeforS6[:-2, 2*(indb-1)+1]), np.absolute(M_BeamsForceS6[:-2, 2*(indb-1)+1]), 'r.-')
        ax.set_xlabel('Curvature (1/m)')
        ax.set_ylabel('Moment (kN*m)')
        ax.set_title('Column')
        ax.grid(True)
        self.ui.MomFi.canvas.draw()
        self.ui.MomFi.canvas.show()

    def IDA(self):
        exec(open("IDA.py").read())


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = MyForm()
    w.show()
    sys.exit(app.exec_())
