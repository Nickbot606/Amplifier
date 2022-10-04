import math
import matplotlib.pyplot as plt

class CSAmp(object):

    def __init__(self,param,debug=False):
        #List key: [vi,Ri,R1,R2,R3,RD,Rs,RL,C1,C2,C3,Vplus,Vmin,K0,Vto,Lambda]
        self.vi = param[0]
        self.Ri = param[1]
        self.R1 = param[2]
        self.R2 = param[3]
        self.R3 = param[4]
        self.RD = param[5]
        self.Rs = param[6]
        self.RL = param[7]
        self.C1 = param[8]
        self.C2 = param[9]
        self.C3 = param[10]
        self.Vplus = param[11]
        self.Vmin = param[12]
        self.K0 = param[13]
        self.Vto = param[14]
        self.Lambda = param[15]

        #Solutions
        self.Sols = self.shortWay(debug)

    def shortWay(self,debug=False):
        Vgg = (self.Vplus*self.R2+(self.Vmin*self.R1))/(self.R1+self.R2)
        Rbb = (self.R1*self.R2)/(self.R1+self.R2)
        V1 = Vgg - self.Vmin - self.Vto
        Id = ((math.sqrt(1 + (4 * self.K0 * V1 * self.Rs)) - 1) / (2 * math.sqrt(self.K0) * self.Rs)) ** 2
        Vds = (self.Vplus - (Id * self.RD) - (self.Vmin + (Id * self.Rs)))
        K = self.K0 * (1 + (self.Lambda * Vds))
        gm = 2 * (math.sqrt(K * Id))
        rs = 1/gm
        r0 = ((1/self.Lambda)+Vds)/(Id)
        rtd = (self.RD * self.RL) / (self.RD + self.RL)
        Vtg = (self.vi*(Rbb))/(Rbb+self.Ri)
        rts = (self.Rs * self.R3) / (self.Rs + self.R3)
        Gmg = (1 / (rs + ((rts * r0)) / (rts + r0))) * (r0 / (r0 + rts))
        rid = r0*(1+(rts/self.Rs))+rts
        vout = -1*Gmg*((rid*rtd)/(rid+rtd))*Vtg
        gain = vout/self.vi
        a = [Vgg,Rbb,V1,Id,K,gm,rs,Vds,r0,rtd,Vtg,rts,Gmg,rid,vout,gain]
        labels = ["Vgg","Rbb","V1","Id","K","gm","rs","Vds","r0","rtd","Vtg","rts","Gmg","rid","vout","gain"]
        if (debug):
            for i in range(0,len(labels)):
                print(labels[i]+" :"+str(a[i]))
        return a

    def getGain(self):
        return self.Sols[15]


class CDAmp(object):
    #Common drain amplifier
    def __init__(self,par,debug=False):
        self.debug = debug
        self.Ri = par[0]
        self.R1 = par[1]
        self.R2 = par[2]
        self.Rd = par[3]
        self.Rs = par[4]
        self.R3 = par[5]
        self.RL = par[6]
        self.Vplus = par[7]
        self.Vmin = par[8]
        self.K0 = par[9]
        self.Vto = par[10]
        self.Lambda = par[11]
        self.Vi = par[12]
        self.Sols = self.smallSignal()

    def smallSignal(self):
        #Solve for Thevenin
        Vgg = (self.Vplus*self.R2+(self.Vmin*self.R1))/(self.R1+self.R2)
        Rbb = (self.R1*self.R2)/(self.R1+self.R2)
        V1 = Vgg - self.Vmin - self.Vto
        Id = ((math.sqrt(1+(4*self.K0*V1*self.Rs))-1)/(2*math.sqrt(self.K0)*self.Rs))**2
        Vds = self.Vplus - (self.Vmin + (Id * self.Rs))
        K = self.K0*(1+(self.Lambda*Vds))
        gm = 2*math.sqrt(K*Id)
        rs = 1/gm
        r0 = ((1/self.Lambda)+(Vds))/(Id)
        vtg = (self.Vi*(Rbb))/(self.Ri+(Rbb))
        rtg = ((self.Ri ** -1) + (self.R1 ** -1) + (self.R2 ** -1)) ** -1
        vs = vtg*(r0/(rs+r0))
        ris = (rs*r0)/(rs+r0)
        v0 = (vs*((self.Rs*self.RL)/(self.Rs+self.RL)))/(ris+((self.Rs*self.RL)/(self.Rs+self.RL)))
        gain = v0/self.Vi
        rin = Rbb
        rout = (ris*self.Rs)/(self.Rs+ris)

        if (self.debug):
            print("\nResults:\n")
            print("Vgg :"+str(Vgg))
            print("Rbb :"+str(Rbb))
            print("Id :"+str(Id))
            print("Vds :"+str(Vds))
            print("gm :"+str(gm))
            print("rs :"+str(rs))
            print("r0 :"+str(r0))
            print("vtg :"+str(vtg))
            print("rtg :"+str(rtg))
            print("vs :"+str(vs))
            print("ris :"+str(ris))
            print("v0 :"+str(v0))
            print("gain! :"+str(gain))
            print("rin :"+str(rin))
            print("rout :"+str(rout))

        return [Vgg,Rbb,Id,Vds,gm,rs,r0,vtg,rtg,vs,ris,v0,gain,rin,rout]

    def getGain(self):
        return self.Sols[12]


class Brute(object):
    #This is to iterate through all the values of a certain aspect
    #iter = iterations, ran = range of the stuff, par = parameters, ind = index of what to replace
    def __init__(self,xlab,ylab,iter,ran,par,ind):
        self.gph = self.populateCS(iter,ran,par,ind,"CSAmp")
        self.xlab = xlab
        self.ylab = ylab

        plt.plot(self.gph[0],self.gph[1])
        plt.xlabel(self.xlab)
        plt.ylabel(self.ylab)
        plt.title(self.xlab+" vs. "+self.ylab)
        #Good morning sunshine, the earth says Hello!
        plt.show()

    def populateCS(self,iter,ran,par,ind,obj):
        a = [[],[]]
        i = ran[0]
        if (obj == "CSAmp"):
            while (i <= ran[1]):
                a[0].append(i)
                par[ind] = i
                i += iter
                a[1].append(CSAmp(par).getGain())
        elif (obj == "CDAmp"):
            while (i <= ran[1]):
                a[0].append(i)
                par[ind] = i
                i += iter
                a[1].append(CDAmp(par).getGain())
        return a


if __name__ == '__main__':
    print("Hello welcome to the amp circuit")
    whichvar = 2
    pars = [1,10,1000000,5000000,0,10000,3000,2000000,1,1,1,3,-3,.001,.5,0.06]
    vardict = ["vin","Ri","R1","R2","R3","RD","Rs","RL","C1","C2","C3","Vplus","Vmin","K0","Vto","Lambda"]
    par = [1,5000,5000000,1000000,50,10000,3000,20000,1,1,1,24,-24,0.001,1.75,0.016]
    #D = CSAmp(pars,True)
    Brute(vardict[whichvar],"Gain",1000,[100,10000000],pars,whichvar)
