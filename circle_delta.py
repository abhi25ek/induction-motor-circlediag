import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import math
import numpy as np
from sympy import *
from sympy.geometry import *

# Prompts the user to enter the details of the 3-ph induction motor
print("Enter the details of the 3-phase induction motor\n")
v = float(input("Enter rated voltage\n"))
p = float(input("Enter number of poles\n"))
freq = float(input("Enter operating frequency\n"))
po = float(input("Enter shaft output power in Watts\n\n"))
srr = float(input("Enter stator to rotor losses ratio\n"))

# Enter no-load data
print("Enter the details of No-load test\n")
vnp = float(input("Enter No-load line voltage in Volts\n"))
inp = float(input("Enter No-load line current in Amperes\n"))
pnp = float(input("Enter No-load power in Watts\n"))

# Enter short-circuit data
print("Enter the details of Blocked test\n")
vbr = float(input("Enter Blocked-rotor line voltage in Volts\n"))
ibr = float(input("Enter Blocked-rotor current in Amperes\n"))
pbr = float(input("Enter Blocked-rotor power in Watts\n"))

# calculations to choose scale of the diagram
tnp = math.acos(pnp/(math.sqrt(3)*vnp*inp))
tbr = math.acos(pbr/(math.sqrt(3)*vbr*ibr))
c_sc = (v*ibr)/vbr


print("\nThe short circuit current is: ")
print(c_sc)


scale = float(
    input("\nChoose an appropriate scale in Amperes per centimetres: "))


def minimum(a, b):

    if a <= b:
        return a
    else:
        return b


class Star:
    def __init__(self, rated_voltage, poles, frequency, power, voltage_np, current_np, power_np, voltage_br, current_br, power_br, stator_rotor_ratio, scale):

        self.text_kwargs = dict(ha='right', va='bottom',
                                fontsize=14, color='#000000')
        self.text_kwargs1 = dict(ha='left', va='bottom',
                                 fontsize=14, color='#000000')

        self.text_kwargs2 = dict(ha='left', va='bottom',
                                 fontsize=14, color='#000000')

        self.rated_voltage = rated_voltage
        self.poles = poles
        self.frequecy = frequency
        self.power = power
        self.voltage_np = voltage_np
        self.voltage_br = voltage_br
        self.current_np = current_np
        self.power_np = power_np
        self.power_br = power_br
        self.current_br = current_br
        self.stator_rotor_ratio = stator_rotor_ratio
        self.scale = scale

        # calculation of no-load and short circuit current per phase and power factors
        self.current_oc = (current_np/self.scale)/math.sqrt(3)

        self.theta_np = (math.acos((self.power_np/3) /
                                   (self.voltage_np*(self.current_np/math.sqrt(3)))))
        self.theta_br = (math.acos((self.power_br/3) /
                                   (self.voltage_br*(self.current_br/math.sqrt(3)))))
        self.current_sc = (
            (self.rated_voltage*(self.current_br/math.sqrt(3)))/self.voltage_br)/self.scale

        self.qnpy = self.current_oc*(math.cos(self.theta_np))
        self.qnpx = self.current_oc*(math.sin(self.theta_np))
        self.dscy = self.current_sc*(math.cos(self.theta_br))
        self.dscx = self.current_sc*(math.sin(self.theta_br))

        # slope of output line
        self.olhorizon = math.atan(
            (self.dscy - self.qnpy)/(self.dscx - self.qnpx))

        self.outmidx = (self.qnpx + self.dscx)/2
        self.outmidy = (self.qnpy + self.dscy)/2

        self.outmiddist = math.sqrt(
            (self.qnpx - self.outmidx)**2+(self.qnpy-self.outmidy)**2)

        # center of the circle
        self.centrex = self.qnpx + (self.outmiddist)/(math.cos(self.olhorizon))

        # radius of the circle
        self.rad = (self.outmiddist)/(math.cos(self.olhorizon))

        self.Fx = self.dscx
        self.Fy = self.qnpy

        # Torque line end point
        self.Gx = (self.stator_rotor_ratio*self.dscx + self.Fx) / \
            (self.stator_rotor_ratio + 1)
        self.Gy = (self.stator_rotor_ratio*self.dscy + self.Fy) / \
            (self.stator_rotor_ratio + 1)

        # powerscale
        self.powerscale_pp = (self.rated_voltage)*self.scale
        self.fullpower_pp = self.power/3

        self.Ddashy = self.dscy + self.fullpower_pp/self.powerscale_pp

        self.c = (self.Ddashy - (math.tan(self.olhorizon)*self.dscx))

        # calculations to the operating point
        self.c1 = Circle(Point(sympify(self.centrex, rational=True), sympify(
            self.qnpy, rational=True)), sympify(self.rad, rational=True))
        self.L1 = Line(Point(sympify(self.dscx, rational=True), sympify(
            self.Ddashy, rational=True)), Point(0, sympify(self.c, rational=True)))

        self.intersectionQ = intersection(self.c1, self.L1)
        (self.Q0a, self.Q0b) = self.intersectionQ[0].evalf()
        (self.Q1a, self.Q1b) = self.intersectionQ[1].evalf()
        self.Q0 = math.sqrt(self.Q0a**2 + self.Q0b**2)
        self.Q1 = math.sqrt(self.Q1a**2 + self.Q1b**2)
        self.Qdash = minimum(self.Q0, self.Q1)
        if self.Qdash == self.Q0:
            self.Q = (self.Q0a, self.Q0b)
        else:
            self.Q = (self.Q1a, self.Q1b)

        # Calculations to find slip, max output power, max torque, slip at max torque
        self.outputline = Line(Point(sympify(self.qnpx, rational=True), sympify(
            self.qnpy, rational=True)), Point(sympify(self.dscx, rational=True), sympify(self.dscy, rational=True)))
        self.torqueline = Line(Point(sympify(self.qnpx, rational=True), sympify(
            self.qnpy, rational=True)), Point(sympify(self.Gx, rational=True), sympify(self.Gy, rational=True)))
        self.qbaseline = Line(Point(sympify(self.qnpx, rational=True), sympify(
            self.qnpy, rational=True)), Point(sympify(self.Fx, rational=True), sympify(self.Fy, rational=True)))
        self.inputQline = Line(Point(sympify(self.Q[0], rational=True), sympify(
            self.Q[1], rational=True)), Point(sympify(self.Q[0], rational=True), 0))

        self.intersectionH = intersection(self.inputQline, self.outputline)
        (self.Hx, self.Hy) = self.intersectionH[0].evalf()
        self.intersectionJ = intersection(self.inputQline, self.torqueline)
        (self.Jx, self.Jy) = self.intersectionJ[0].evalf()
        self.intersectionK = intersection(self.inputQline, self.qbaseline)
        (self.Kx, self.Ky) = self.intersectionK[0].evalf()

        self.statorI2R = abs(self.Ky - self.Jy)
        self.rotorI2R = abs(self.Jy - self.Hy)
        self.rated_slip = abs(self.Hy-self.Jy)/(abs(self.Q[1]-self.Jy))
        self.sychronous_speed = (4*math.pi*self.frequecy)/(self.poles)
        self.airgap_rated = abs(self.Q[1]-self.Jy)
        self.total_input_rated = self.Q[1]
        self.torque_rated = (
            self.airgap_rated/self.sychronous_speed)*self.powerscale_pp*3
        self.efficiency_rated = abs(self.Q[1]-self.Hy)/self.Q[1]

        self.maxoutputangle = math.pi - ((math.pi/2) - self.olhorizon)
        self.cmaxoutput = self.qnpy - \
            math.tan(self.maxoutputangle)*self.centrex
        self.maxoutputline = Line(Point(sympify(self.centrex, rational=True), sympify(
            self.qnpy, rational=True)), Point(0, sympify(self.cmaxoutput, rational=True)))

        self.intersectionA = intersection(self.c1, self.maxoutputline)
        (self.Aax, self.Aay) = self.intersectionA[0].evalf()
        (self.Abx, self.Aby) = self.intersectionA[1].evalf()
        if self.Aay >= 0:
            self.A = (self.Aax, self.Aay)
        else:
            self.A = (self.Abx, self.Aby)

        self.maxoutputlineintersect = Line(Point(sympify(self.A[0], rational=True), sympify(
            self.A[1], rational=True)), Point(sympify(self.A[0], rational=True), 0))

        self.intersectionA1 = intersection(
            self.maxoutputlineintersect, self.outputline)
        (self.A10x, self.A11y) = self.intersectionA1[0].evalf()

        self.maxoutput = abs(self.A[1] - self.A11y)

        self.torquelinemidx = (self.qnpx + self.Gx)/2
        self.torquelinemidy = (self.qnpy + self.Gy)/2

        self.torquelineangle = math.atan(
            (self.Gy - self.qnpy)/(self.Gx - self.qnpx))

        self.maxtorqueangle = math.pi - ((math.pi/2) - self.torquelineangle)

        self.cmaxtorque = self.qnpy-math.tan(self.maxtorqueangle)*self.centrex

        self.maxtorqueline = Line(Point(sympify(self.centrex, rational=True), sympify(
            self.qnpy, rational=True)), Point(0, sympify(self.cmaxtorque, rational=True)))

        self.intersectionB = intersection(self.c1, self.maxtorqueline)
        (self.Bax, self.Bay) = self.intersectionB[0].evalf()
        (self.Bbx, self.Bby) = self.intersectionB[1].evalf()
        if self.Bay >= 0:
            self.B = (self.Bax, self.Bay)
        else:
            self.B = (self.Bbx, self.Bby)

        self.maxtorquelineintersect = Line(Point(sympify(self.B[0], rational=True), sympify(
            self.B[1], rational=True)), Point(sympify(self.B[0], rational=True), 0))

        self.intersectionB1 = intersection(
            self.torqueline, self.maxtorquelineintersect)
        (self.B1x, self.B1y) = self.intersectionB1[0].evalf()

        self.intersectionf1 = intersection(
            self.maxtorquelineintersect, self.outputline)
        (self.f1x, self.f1y) = self.intersectionf1[0].evalf()

        self.maxtorque = (abs(self.B[1] - self.B1y) /
                          self.sychronous_speed)*self.powerscale_pp*3
        self.maxtorqueslip = (abs(self.f1y - self.B1y)) / \
            (abs(self.B[1] - self.B1y))

        self.startingtorque = (abs(self.dscy - self.Gy) /
                               self.sychronous_speed)*self.powerscale_pp*3
        self.inputpowerfactor_rated = math.cos(math.atan(self.Q[0]/self.Q[1]))
        self.rotorpowerfactor_rated = math.cos(
            math.tan((self.Q[0] - self.qnpx)/(self.Q[1] - self.qnpy)))
        self.inputstatorcurrent = math.sqrt((self.Q[0])**2 + (self.Q[1]**2))
        self.inputrotorcurrent = math.sqrt(
            (self.Q[0]-self.qnpx)**2 + (self.Q[1]-self.qnpy)**2)

    def draw_star(self):
        basex = np.array([0, self.centrex + self.rad, 0, 0])
        basey = np.array([0, 0, 0, self.rad + self.qnpy])
        ie_x = np.array([0, self.qnpx])
        ie_y = np.array([0, self.qnpy])
        plt.arrow(0, 0, 0, self.rad + self.qnpy, head_width=1/(self.rad + self.qnpy),
                  head_length=2/(self.rad + self.qnpy), fc='#0099ff', ec='#0099ff', color='#0099ff')
        plt.plot(basex, basey, color='#0099ff')
        plt.plot(ie_x, ie_y, color='#A3E4D7')
        isc_x = np.array([0, self.dscx])
        isc_y = np.array([0, self.dscy])
        plt.plot(isc_x, isc_y, color='#A3E4D7')
        outputlx = np.array([self.qnpx, self.dscx])
        outputly = np.array([self.qnpy, self.dscy])
        plt.plot(outputlx, outputly)
        qbasex = np.array([self.qnpx, self.rad+self.centrex])
        qbasey = np.array([self.qnpy, self.qnpy])
        plt.plot(qbasex, qbasey, color='#E74C3C')
        theta = np.linspace(0, np.pi, 1000)
        circlex = self.centrex + self.rad*np.cos(theta)
        circley = self.qnpy + self.rad*np.sin(theta)
        plt.plot(circlex, circley, color='#E74C3C')
        ohmsclossx = np.array([self.Fx, self.dscx])
        ohmsclossy = np.array([self.Fy, self.dscy])
        plt.plot(ohmsclossx, ohmsclossy)
        torquex = np.array([self.qnpx, self.Gx])
        torquey = np.array([self.qnpy, self.Gy])
        plt.plot(torquex, torquey, color='#6495ED')
        dashx = np.array([self.dscx, self.dscx])
        dashy = np.array([self.dscy, self.Ddashy])
        plt.plot(dashx, dashy)
        Qx = np.array([self.Q[0], self.dscx])
        Qy = np.array([self.Q[1], self.Ddashy])
        plt.plot(Qx, Qy, color='#808B96')
        inputQx = np.array([self.Q[0], self.Q[0]])
        inputQy = np.array([self.Q[1], 0])
        plt.plot(inputQx, inputQy, color='#17A589')
        inputcurrentQx = np.array([0, self.Q[0]])
        inputcurrentQy = np.array([0, self.Q[1]])
        plt.plot(inputcurrentQx, inputcurrentQy, color='#8E44AD')
        mox = np.array([self.A[0], self.A10x])
        moy = np.array([self.A[1], self.A11y])
        plt.plot(mox, moy)
        mox1 = np.array([self.centrex, self.A[0]])
        mox2 = np.array([self.qnpy, self.A[1]])
        plt.plot(mox1, mox2)
        mtx = np.array([self.B[0], self.B[0]])
        mty = np.array([self.B[1], self.B1y])
        plt.plot(mtx, mty)
        mt1x = np.array([self.centrex, self.B[0]])
        mt1y = np.array([self.qnpy, self.B[1]])
        plt.plot(mt1x, mt1y)
        Irotor_ratedx = np.array([self.qnpx, self.Q[0]])
        Irotor_ratedy = np.array([self.qnpy, self.Q[1]])
        plt.plot(Irotor_ratedx, Irotor_ratedy)
        plt.axis("equal")
        plt.text(self.qnpx, self.qnpy, 'Q', **self.text_kwargs)
        plt.text(self.Q[0], self.Q[1], 'P', **self.text_kwargs)
        plt.text(self.dscx, self.dscy, 'D', **self.text_kwargs1)
        plt.text(self.Gx, self.Gy, 'G', **self.text_kwargs1)
        plt.text(self.Fx, self.Fy, 'F', **self.text_kwargs1)
        plt.text(self.dscx, self.Ddashy, 'D\'', **self.text_kwargs1)
        plt.text(self.centrex, self.qnpy, 'C', **self.text_kwargs1)
        plt.text(self.A[0], self.A[1], 'A', **self.text_kwargs1)
        plt.text(self.A10x, self.A11y, 'A1', **self.text_kwargs)
        plt.text(self.B[0], self.B[1], 'B', **self.text_kwargs1)
        plt.text(self.B1x, self.B1y, 'B1', **self.text_kwargs1)
        plt.text(self.f1x, self.f1y, 'f', **self.text_kwargs1)
        plt.text(self.Hx, self.Hy, 'H', **self.text_kwargs1)
        plt.text(self.Jx, self.Jy, 'J', **self.text_kwargs1)
        plt.text(self.Kx, self.Ky, 'K', **self.text_kwargs2)
        plt.text(0, self.rad + self.qnpy, 'V', **self.text_kwargs)
        plt.show()


ob = Star(v, p, freq, po, vnp, inp, pnp, vbr, ibr, pbr, srr, scale)
print("\nNo-load power factor angle is: ")
print((180*ob.theta_np)/math.pi)
print("\nShort-circuit power factor angle is: ")
print((180*ob.theta_br)/math.pi)
print("\nNo-load current is: ")
print(ob.current_np)
print("\nShort circuit current at rated stator voltage is: ")
print(ob.current_sc*ob.scale)
print("\nRadius of the circle is: ")
print(ob.rad)
print("\nLength of DD' is: ")
print(ob.Ddashy - ob.dscy)
print("\nInput stator current at rated conditions is: ")
print(ob.inputstatorcurrent*ob.scale)
print("\nInput rotor current at rated conditions is: ")
print(ob.inputrotorcurrent*ob.scale)
print("\nInput stator pf at rated conditions is: ")
print(ob.inputpowerfactor_rated)
print("\nInput rotor pf at rated conditions is: ")
print(ob.rotorpowerfactor_rated)
print("\nEfficiency at rated conditions is: ")
print(ob.efficiency_rated)
print("\nSlip at rated conditions is: ")
print(ob.rated_slip)
print("\nSynchronous speed is: ")
print(ob.sychronous_speed)
print("\nStarting torque is: ")
print(ob.startingtorque)
print("\nTorque at rated conditions is: ")
print(ob.torque_rated)
print("\nAirgap power at rated conditions is: ")
print(ob.airgap_rated*ob.powerscale_pp*3)
print("\nMaximum Torque is: ")
print(ob.maxtorque)
print("\nSlip at max torque is: ")
print(ob.maxtorqueslip)
print("\nMaximum Output is: ")
print(ob.maxoutput*ob.powerscale_pp*3)
ob.draw_star()
