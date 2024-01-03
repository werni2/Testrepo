import re
import math
from mymath import *
import matplotlib.pyplot as plt
import numpy as np
import random

class AccSection:

    def __init__(self, vi, vf, ai, af, A, D, J):

        self.__s0	= 0.0
        self.__v0	= vi  
        self.__v3b	= vf
        self.__a0	= ai
        self.__a3b	= af

        if 	IsEqual(vi, vf) and IsEqual(ai, af): 

            self.__Acenter = ai
            
            self.__a1a = ai
            self.__v1a = vi
            self.__s1a = self.__s0
            
            self.__a3a = af
            self.__v3a = vf
            
            self.__isTrapez = False
            self.__tau0a = 0.0
            self.__tau2b = 0.0
            self.__tau1  = 0.0
            self.__tau0b = 0.0
            self.__tau2a = 0.0
            
        elif ai <= A and ai >= -D and af <= A and af >= -D: 

            delta_a = af - ai
            delta_v = vf - vi
            if delta_a > 0.0: 
            
                if ai > 0.0: 
                    delta_vm = (af * af - ai * ai) / (2.0 * J)
                elif af < 0.0: 
                    delta_vm = (af * af - ai * ai) / (2.0 * J)
                else:
                    delta_vm =(af * af) / (2.0 * J) - (ai * ai) / (2.0 * J)
                
            else:

                if af > 0.0: 
                    delta_vm = (ai * ai - af * af) / (2.0 * J)
                elif ai < 0.0:    
                    delta_vm = (ai * ai - af * af) / (2.0 * J)
                else:
                    delta_vm = (ai * ai) / (2.0 * J) - (af * af) / (2.0 * J)
                      
            # acceleration?
            self.__dir = 1.0 if delta_v - delta_vm > 0.0 else -1.0
            if self.__dir > 0.0: 
                # yes!
                self.__Acenter = A
                self.__J = J
            else:
                # no!
                self.__Acenter = -D
                self.__J = -J    

            self.__hasPrePhase = ai * self.__dir < 0.0
            if self.__hasPrePhase: 
                # hasPrePhase
                self.__tau0a = -ai / self.__J
                delta_v = 0.5 * ai * self.__tau0a
                delta_s = (vi - (1.0 / 3.0) * (ai * ai) / self.__J) * self.__tau0a
                self.__a1a = 0.0
                self.__v1a = vi + delta_v
                self.__s1a = self.__s0 + delta_s
            else:
                self.__tau0a = 0.0
                self.__a1a = ai
                self.__v1a = vi
                self.__s1a = self.__s0
            
            self.__hasPostPhase = af * self.__dir < 0.0
            if self.__hasPostPhase: 
            
                # hasPostPhase
                self.__tau2b = -self.__a3b / self.__J
                delta_v = 0.5 * self.__a3b * self.__tau2b
                self.__a3a = 0.0
                self.__v3a = self.__v3b - delta_v
            else:
                self.__tau2b = 0.0
                self.__a3a = af
                self.__v3a = vf
            
            if self.__dir * (1.0 if self.__v3a - self.__v1a > 0.0 else -1.0) > 0.0: 
                tmp0 = 1.0 / (2.0 * self.__J)
                tmp1 = 1.0 / (2.0 * self.__J)
                Z = self.__v3a - self.__v1a + tmp1 * (self.__a1a * self.__a1a) + tmp0 * (self.__a3a * self.__a3a)
                N = tmp0 + tmp1
                F = Z / N
                F = Threshold(F, 1E-12)
                if F >= 0.0:
                    F = self.__dir * math.sqrt(F)
                    self.__isTrapez = self.__dir * F > self.__dir * self.__Acenter
                    if self.__isTrapez:
                        self.__a1b = self.__a2 = self.__Acenter
                    else:
                        self.__Acenter = self.__a1b = self.__a2 = F
                        self.__tau1 = 0.0
                    
                    self.__tau0b = (self.__a1b - self.__a1a) / self.__J 
                    self.__tau2a = (self.__a2 - self.__a3a) / self.__J 
                else:
                    raise Exception("Acceleration must not be negative")

            self.__v1b 	= self.__v1a + (self.__a1a + 0.5 * self.__J * self.__tau0b) * self.__tau0b 
            self.__v2 	= self.__v3a - (self.__a2 - 0.5 * self.__J * self.__tau2a) * self.__tau2a if self.__isTrapez else self.__v1b 

            if self.__isTrapez:
                self.__tau1 = (self.__v2 - self.__v1b) / self.__a2

            # compute pos bounderies
            self.__s1b 	= self.__s1a + (self.__v1a + (0.5 * self.__a1a + (1.0 / 6.0) * self.__J * self.__tau0b) * self.__tau0b) * self.__tau0b
            self.__s2 	= self.__s1b + (self.__v1b + 0.5 * self.__a1b * self.__tau1) * self.__tau1 if self.__isTrapez else self.__s1b
            self.__s3a 	= self.__s2 + (self.__v2 + (0.5 * self.__a2 - (1.0 / 6.0) * self.__J * self.__tau2a) * self.__tau2a) * self.__tau2a
            self.__s3b 	= self.__s3a + (self.__v3a + (0.5 * self.__a3a - (1.0 / 6.0) * self.__J * self.__tau2b) * self.__tau2b) * self.__tau2b

            # compute time bounderies
            self.__t0 	= 0.0
            self.__t1a 	= self.__t0 + self.__tau0a 
            self.__t1b 	= self.__t1a + self.__tau0b 
            self.__t2 	= self.__t1b + self.__tau1 
            self.__t3a 	= self.__t2 + self.__tau2a 
            self.__t3b 	= self.__t3a + self.__tau2b 

            # set base data
            self.Duration   = self.__t3b
            self.Length     = self.__s3b

    def GetMotionState(self, t):

        if IsLowerOrEqual(t, self.__t0):
            pos =  0.0
            vel =  self.__v0
            acc =  self.__a0
        elif self.__hasPrePhase and IsLowerOrEqual(t, self.__t1a):
            # phase 0
            t  = t - self.__t0
            pos =  self.__s0 + (self.__v0 + (0.5 * self.__a0 + (1.0 / 6.0) * self.__J * t) * t) * t
            vel =  self.__v0 + (self.__a0 + 0.5 * self.__J * t) * t
            acc =  self.__a0 + self.__J * t
        elif IsLowerOrEqual(t, self.__t1b):
            # phase 0
            t = t - self.__t1a
            pos =  self.__s1a + (self.__v1a + (0.5 * self.__a1a + (1.0 / 6.0) * self.__J * t) * t) * t
            vel =  self.__v1a + (self.__a1a + 0.5 * self.__J * t) * t
            acc =  self.__a1a + self.__J * t
        elif self.__isTrapez and IsLowerOrEqual(t, self.__t2):
            # phase 1
            t = t - self.__t1b
            pos =  self.__s1b + (self.__v1b + 0.5 * self.__Acenter * t) * t 
            vel =  self.__v1b + self.__Acenter * t
            acc =  self.__Acenter
        elif IsLowerOrEqual(t, self.__t3a):
            # phase 2
            t = t - self.__t2
            pos =  self.__s2 + (self.__v2 + (0.5 * self.__Acenter - (1.0 / 6.0) * self.__J * t) * t) * t
            vel =  self.__v2 + (self.__a2 - 0.5 * self.__J * t) * t
            acc =  self.__Acenter - self.__J * t
        elif self.__hasPostPhase and IsLowerOrEqual(t, self.__t3b):
            # phase 0
            t = t - self.__t3a
            pos =  self.__s3a + (self.__v3a + (0.5 * self.__a3a - (1.0 / 6.0) * self.__J * t) * t) * t
            vel =  self.__v3a + (self.__a3a - 0.5 * self.__J * t) * t
            acc =  self.__a3a - self.__J * t
        else:
            pos =  self.__s3b
            vel =  self.__v3b
            acc =  self.__a3b

        return pos, vel, acc
            
class MotionProfile:

    def __init__(self, vi, vf, ai, af, V, A, D, J, L):

        def __bisecFunc1(V):
            self.__accSec = AccSection(vi, V, ai, 0, A, D, J)
            self.__decSec = AccSection(V, vf, 0, af, A, D, J)
            return L - self.__accSec.Length - self.__decSec.Length
        
        def __bisecFunc2(ac):
            self.__accSec = AccSection(vi, 0, ai, ac, A, D, J)
            self.__decSec = AccSection(0, vf, ac, af, A, D, J)
            return L - self.__accSec.Length - self.__decSec.Length
        
        def __bisecFunc3(vsplit):
            self.__accSec = AccSection(vi, vsplit, ai, A, A, D, J)
            self.__decSec = AccSection(vsplit, vf, A, af, A, D, J)
            return L - self.__accSec.Length - self.__decSec.Length
        
        def __bisecFunc4(vsplit):
            self.__accSec = AccSection(vi, vsplit, ai, -D, A, D, J)
            self.__decSec = AccSection(vsplit, vf, -D, af, A, D, J)
            return L - self.__accSec.Length - self.__decSec.Length

        vi = Threshold(vi, 1E-12)
        vf = Threshold(vf, 1E-12)
        ai = Threshold(ai, 1E-12)
        af = Threshold(af, 1E-12) 

        if vi + 0.5 * ai * abs(ai) / J < 0.0:
            swap(A, D)

        self.__accSec = AccSection(vi, vf, ai, af, A, D, J) 
        if IsEqual(self.__accSec.Length, L):
            self.__decSec = AccSection(vf, vf, af, af, A, D, J)
        elif L > self.__accSec.Length:
            # normal path
            diff = __bisecFunc1(V)
            if IsGreaterOrEqual(diff, 0.0):
                self.__Vcenter = V
                self.__TMitte = Threshold(diff / self.__Vcenter, 1E-12)
            elif L * self.__accSec.Length > 0:
                self.__Vcenter = bisection(0, V, __bisecFunc1)
                self.__TMitte = 0.0
            else:
                fA = __bisecFunc2(-D)
                fB = __bisecFunc2(A)
                if fA * fB < 0.0:
                    bisection(-D, A, __bisecFunc2)
                    self.__Vcenter = 0
                    self.__TMitte = 0.0
                else:
                    bisection(-V, V, __bisecFunc3)
                    self.__Vcenter = 0
                    self.__TMitte = 0.0
        else:
            # reversing
            diff = __bisecFunc1(-V)
            if IsLowerOrEqual(diff, 0.0):
                self.__Vcenter = -V
                self.__TMitte = Threshold(diff / self.__Vcenter, 1E-12)
            elif L * self.__accSec.Length > 0:
                self.__Vcenter = bisection(-V, 0, __bisecFunc1)
                self.__TMitte = 0.0
            else:
                fA = __bisecFunc2(-D)
                fB = __bisecFunc2(A)
                if fA * fB < 0.0:
                    bisection(-D, A, __bisecFunc2)
                    self.__Vcenter = 0
                    self.__TMitte = 0.0
                else:
                    bisection(-V, V, __bisecFunc4)
                    self.__Vcenter = 0
                    self.__TMitte = 0.0

        self.xValid = True
        self.Duration = self.__accSec.Duration + self.__decSec.Duration + self.__TMitte
        self.Length = L  

    def GetMotionState(self, t):

        if t < self.__accSec.Duration:
            pos, vel, acc = self.__accSec.GetMotionState(t)
        elif t < self.__accSec.Duration + self.__TMitte:
            t = t - self.__accSec.Duration
            pos = self.__accSec.Length + self.__Vcenter * t
            vel = self.__Vcenter
            acc = 0.0
        elif IsLowerOrEqual(t, self.__accSec.Duration + self.__TMitte + self.__decSec.Duration):
            t = t - self.__accSec.Duration - self.__TMitte
            pos, vel, acc = self.__decSec.GetMotionState(t)
            pos = pos + self.__accSec.Length + self.__Vcenter * self.__TMitte 
        else:
            raise Exception("time out of range")
        
        return pos, vel, acc

# Test
if __name__ == "__main__":
    print("main")

    #profile = AccSection(100, -100, 100, 2100, 5000, 5000, 10000)
    #profile = MotionProfile(100, -100, 100, 200, 5000, 5000, 5000, 10000, 0) 
    while True:

        V = random.uniform(100, 5000)
        A = random.uniform(1000, 5000)
        D = random.uniform(1000, 5000)
        J = random.uniform(1, 1e10)
        L = random.uniform(-10000, 10000)
        vi = random.uniform(-V, V)
        vf = random.uniform(-V, V)
        ai = random.uniform(-D, A)
        af = random.uniform(-D, A)

        profile = MotionProfile(vi, vf, ai, af, V, A, D, J, L) 

        sfinal, vfinal, afinal = profile.GetMotionState(profile.Duration) 
        if not (IsEqual(L, sfinal, 1e-6) and IsEqual(vf, vfinal, 1e-6) and IsEqual(af, afinal, 1e-6)):
            raise Exception("invalif final values")
            

        t = []
        s = []
        v = []
        a = []

        for time in np.linspace(0, profile.Duration, 1000):
            pos, vel, acc = profile.GetMotionState(time) 
            t.append(time)
            s.append(pos)
            v.append(vel)
            a.append(acc)
        
        # Erstellen Sie Subplots (3 Zeilen, 1 Spalte)
        plt.subplot(3, 1, 1)
        plt.plot(t, s)
        plt.title('distance')

        plt.subplot(3, 1, 2)
        plt.plot(t, v)
        plt.title('velocity')

        plt.subplot(3, 1, 3)
        plt.plot(t, a)
        plt.title('acceleration')

        # Einstellungen für das gesamte Diagramm
        plt.xlabel('time')
        plt.tight_layout()  # Verbessert den Layout-Abstand zwischen Subplots

        plt.show()

        #print("")

        # blablabla

    print("")