from time import sleep
from math import ceil
from os import path
import numpy as np
import textwrap
import operator

def findfreq(freq,Ecc,maxsamp,factor):
    #the single-tone algorithm uses a repeating waveform to reduce the length of waveforms needed to program the AWG.
    #This repeating waveform discretizes the allowed frequency spectrum of the output.
    #Therefore the number of oscillations per repeating segment, and the duration of this segment need to be chosen carefully in order to achieve an accurate frequency.
    #This function accomplishes this using 2 criteria. There is a maxsamp, which quantifies the maximum number of samples used for each segment.
    #And Ecc, which quantifies the desired accuracy of the achieved frequency.
    #If the desired accuracy can be reached with less than the maximum number of samples. The minimum number of samples which satisfy the desired accuracy will be chosen.
    #If this first criterion cannot be achieved. The algorithm will choose the parameters for the frequency closest to the desired frequency achievable within the maxsamp.
    q=16*freq/(24e8/2**factor)
    qaprox=1e100
    a,b,c,d=0,1,1,1
    e,f=1,1
    i=1
    while abs((24e8/(2**factor))*qaprox)/16>freq/Ecc and f*16<maxsamp:
        e,f=a+c,b+d
        if f*16<maxsamp:
            if q>e/f:
                a,b=e,f
            else:
                c,d=e,f
        qaprox=min(abs((q-a/b)),abs((q-c/d)))
        
        
        i+=1
    if abs((q-a/b))<abs((q-c/d)):
        e,f=a,b
    else:
        e,f=c,d
    factor2=int(100/(f*16))+1
    WN1=e*factor2
    LEN1=f*16*factor2
    

    freqout=WN1/LEN1*(24e8/2**factor)
    
    return(WN1,LEN1,freqout)
'''
def findfreq(freq,Ecc,maxsamp,factor):
    #the single-tone algorithm uses a repeating waveform to reduce the length of waveforms needed to program the AWG.
    #This repeating waveform discretizes the allowed frequency spectrum of the output.
    #Therefore the number of oscillations per repeating segment, and the duration of this segment need to be chosen carefully in order to achieve an accurate frequency.
    #This function accomplishes this using 2 criteria. There is a maxsamp, which quantifies the maximum number of samples used for each segment.
    #And Ecc, which quantifies the desired accuracy of the achieved frequency.
    #If the desired accuracy can be reached with less than the maximum number of segments. The minimum number of segments which satisfy the desired accuracy will be chosen.
    #If this first criterion cannot be achieved. The algorithm will choose the parameters for the frequency closest to the desired frequency achievable within the maxsamp.
    tSamp=1/(2.4e9)*2**factor
    SampPerPeriod=1/(freq*tSamp)
    maxPeriods=int(((maxsamp/SampPerPeriod)))
    

    #Create a interger numpy-array counting up to the maximum number of loops.
    Ps=np.linspace(1,maxPeriods,maxPeriods,dtype=int)

    #First check which must be satisfied is that the segment contains a multiple of 16 samples. This is because of HDAWG playback criteria.
    #If this check is not met, then the segments would be zero-extended to a multiple of 16 samples. 
    #This check also ensures the repeating segments are longer than 100 samples, as the AWG has a minimum waveform length.
    Lfloor=((Ps/(freq*tSamp))%16<0.99)*(Ps/(freq*tSamp)>100)
    Lceil=((Ps/(freq*tSamp))%16>15.01)*(Ps/(freq*tSamp)>100)
    if True not in Lc and True not in Lfloor: #if these requirements cannot be met
        raise(Exception("increase the maxsamp"))

    #The second check validates whether each loop number has an achievable frequency close enough to the desired frequency.
    remainder=np.divide((Ps/freq%(tSamp)),Ps)
    Kfloor=remainder<1/(freq*Ecc)
    Kceil=np.abs(np.divide(tSamp,Ps)-remainder)<1/(freq*Ecc)

    #each of these checks produce a truth array, (False,False,True...etc.). Which is multiplied to form two overall truth tables,
    #depending on whether the desired frequency is under the number of samples calculated(ceil) or over the number of samples calculated(floor)
    Kfloor=Lfloor*Kfloor
    Kceil=Lceil*Kceil

    #A&B verify if the first criterion can be reached. 
    A,B= True in Kfloor, True in Kceil
    
    if A and B: #if both ceil and floor can reach first criterion, which can reach it faster?
        Kfloor=np.where(Kfloor)
        WN1f=Kfloor[0][0]+1
        Kceil=np.where(Kceil)
        WN1c=Kceil[0][0]+1
        if WN1c<WN1f:
            C='c'
            WN1=WN1c
        else:
            C='f'
            WN1=WN1f

    #if only one works, take what works.
    elif A:
        K=np.where(Kfloor)
        WN1=K[0][0]+1
        C='f'

    elif B:
        K=np.where(Kceil)
        WN1=K[0][0]+1
        C='c'
    else:#else, use the second criterion.
        Pf=remainder
        
        Pc=np.abs(np.divide(tSamp,Ps)-remainder)
        A=np.min(Pf+np.logical_not(Lfloor)*10000) 
        B=np.min(Pc+np.logical_not(Lceil)*10000)
        if B>A:
            WN1=np.where(Pf==np.min(Pf+np.logical_not(Lfloor)*1000))[0][0]+1
            C='f'
        else:
            WN1=np.where(Pc==np.min(Pc+np.logical_not(Lceil)*1000))[0][0]+1
            C='c'
    #Find the ultimate lenght of this segment, depending on whether the flooring or ceiling criterion is used.
    if C=='f':
        LEN1=int(WN1*2.4e9/(freq*2**factor))
    else:
        LEN1=ceil(WN1*2.4e9/(freq*2**factor))
    
    freqout=WN1/LEN1*24e8/(2**factor)
    
    if abs(freqout-freq)>freq/Ecc:
        print("Sufficiently accurate frequency could not be found, most accurate frequency available is off by "+ str(freqout-freq) )
    else:
        print("A sufficiently accurate frequency could be found, its value is "+str(freqout)+". the desired frequency was "+str(freq))
    return WN1,LEN1,freqout
'''

def HDAWGcalibration(amp1=0.5,amp2=0.5,freq=20e6,phase1=0,phase2=np.pi/2):
        damp,dphase=0,0
        #Changes the nominal amplitude and phase of channel 2 in order to improve IQ modulation using the R&S device. Calibrated to the R&S Device.
        if amp2==amp1:
            damp=amp1*amp1*0.0021+0.0044*amp1+0.0003
        if phase2==phase1+np.pi/2 or phase2==phase1-np.pi/2:
            dphase=-freq*1e-6*0.0008+0.0006
        return damp,dphase

class Sequence():
    #this class defines and builds the sequences which are loaded onto the HDAWG.
    def __init__(self,Trigchan:int=1,While:bool=True, Ecc:float=1e7, maxsamp:int=1e7, acc:float=50, t0=0):
        #This initializes the sequence, and defines a few global parameters used for defining the accuracy of the tones, and defines whether the sequence should loop. ("While")
        
        self.seq = textwrap.dedent(
            """
            while(1){
            waitDigTrigger(DTR);
            }
            """
        )
        
        if While!=True:
            self.seq = self.seq.replace("while(1){","").replace("}","")
            self.While=False
        else:
            self.While=True
        #defines which channel should trigger the sequence.
        if Trigchan==False:
            self.seq = self.seq.replace("waitDigTrigger(DTR);", "")
        self.seq = self.seq.replace("DTR", str(Trigchan)) 

        #defines various values related to the desired frequency accuracy of the HDAWG output.
        self.Ecc = Ecc
        self.msamp  = maxsamp
        self.acc = acc
        self.t=t0
        #prepares a dictionary which will store the approximated frequencies for the looping segments computed by freqout(). This serves to reduce redundancy.
        self.freqdict={}


    def modulator_on(self):
        #Test code, should turn on modulation on channel 1 and 2.
        #It seems to fail.
        if self.While==True:
            self.seq= self.seq[:len(self.seq)-2]

        self.seq+=textwrap.dedent(
            """
            setInt("modulation/0/modulation_mode(0)", 1);
            setInt("modulation/0/modulation_mode(1)", 2);
            
            """
        )
        if self.While==True:
            self.seq +="}"
        
    
    def modulator_off(self):
        #Test code, should turn off modulation on channel 1 and 2.
        #It seems to fail.
        if self.While==True:
            self.seq= self.seq[:len(self.seq)-2]

        self.seq+=textwrap.dedent(
            """
            setInt("modulation/0/modulation_mode(0)", 0);
            setInt("modulation/0/modulation_mode(1)", 0);
            
            """
        )
        if self.While==True:
            self.seq +="}"

    

    
    def addsingletone( #This function produces a single-tone two-channel signal.
        self,
        dur:float       = 0.001,
        freq:float      = 1e6,
        amp:float       = 1.0,
        amp2:float      = None,
        phase1:float    = None,
        phase2:float    = None,
        Mark:list       = [0,0,0,0],#if and where a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto phase1 and phase2.
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2  #2nd Output channel
    ):
        Flag1,Flag2=False,False
        if amp2==None:
            amp2=amp
            Flag1=True

        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
            Flag2=True
        self.t+=dur

        #default phases of channel 1 and 2 are sine and cosine.
        if phase1==None:
            phase1=0
        if phase2==None:
            phase2=phase1-np.pi/2

        damp,dphase=HDAWGcalibration(amp,amp2,freq,phase1,phase2)
        if Flag1:
            amp2+=damp
        if Flag2:
            phase2+=dphase
        
        #For various calcs, need the following constants.
        per=freq*dur
        if dur>1e-7:
            factor=max(int(np.log2(2.4e9/(self.acc*freq))),0)
        else:
            factor=0
        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        if freq in self.freqdict:
            [WN1,LEN1,freqout,factor2]=self.freqdict[freq]
            if dur>1e-7 and factor!=factor2:
                WN1,LEN1,freqout=findfreq(freq,self.Ecc,self.msamp,factor)
                factor=factor2
        else:
            WN1,LEN1,freqout=findfreq(freq,self.Ecc,self.msamp,factor)
            self.freqdict[freq]=[WN1,LEN1,freqout,factor]
        N=int(round(2.4e9*dur*1/(2**factor),5))

        #Interpreting part of the output of freqout()
        if WN1==0:
            per=0
        else:
            per=int(per/WN1)
        
        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addsingletone does not)
        phase1+=(freqout*2*np.pi*phasetime)%(2*np.pi)
        phase2+=(freqout*2*np.pi*phasetime)%(2*np.pi)

        #Find the length of the final, unrepeated, waveform. Checking to ensure the length of this waveform is positive.
        LEN2=N-LEN1*per
        while(LEN2<0):
            per=per-1
            LEN2=N-LEN1*per
        per2=freqout*(LEN2/24e8)*2**factor
        
        #Makes sure the final waveform is long enough
        if round(LEN2,5)<100 and per>0:
            per-=1
            LEN2+=LEN1
            per2+=WN1

        factor2=factor
        while LEN2%16!=0 and factor2>0:
            factor2-=1
            LEN2*=2
        if LEN2%16!=0:
            print("The precision required by this pulse is not achievable. Will zero extend this part of the pulse, which may cause various issues.")



        #Starting the awg_program for this section depending on whether there should be a repeating section:
        if per==0:
            awg_program = textwrap.dedent(
                """
                    playWave(1,marker(LEN2,1)+WAFE1,2,marker(LEN2,1)+WAFE2,3,marker(LEN2,1)+WAFE3,4,marker(LEN2,1)+WAFE4,FACT2);
                """
            )
        else:
            awg_program = textwrap.dedent(
                """
                    repeat(PER){
                        playWave(1,marker(LEN1,1)+WAVE1,2,marker(LEN1,1)+WAVE2,3,marker(LEN1,1)+WAVE3,4,marker(LEN1,1)+WAVE4,FACT1);
                    }
                    playWave(1,marker(LEN2,1)+WAFE1,2,marker(LEN2,1)+WAFE2,3,marker(LEN2,1)+WAFE3,4,marker(LEN2,1)+WAFE4,FACT2);
                """
            )

        
        

        #in order to have arbitrary-channel outputs, the following code makes sure that the outputs: OUT1 & OUT2 end up in the correct spots.
        #This is not needed in the current state of the experiment.

        OUT1="+AMP1*sine(LEN1, P1, WN1)"
        OUT2="+AMP2*sine(LEN1, P2, WN1)"

        
        
        if CHAN1==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT1)
        else:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT2)
        awg_program=awg_program.replace("+WAVE1","")
        awg_program=awg_program.replace("+WAVE2","")
        awg_program=awg_program.replace("+WAVE3","")
        awg_program=awg_program.replace("+WAVE4","")
        

        awg_program=awg_program.replace("WAFE","WAVE")
        OUT1="+AMP1*sine(LEN2, P1,WN2)"
        OUT2="+AMP2*sine(LEN2, P2,WN2)"

        if CHAN1==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT1)
        else:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT2)
        awg_program=awg_program.replace("+WAVE1","")
        awg_program=awg_program.replace("+WAVE2","")
        awg_program=awg_program.replace("+WAVE3","")
        awg_program=awg_program.replace("+WAVE4","")
        
        


        
        #The following code makes sure that the correct markers are retained, while removing all marker waveforms which are unneeded.        
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1)+",str(i+1)+",")




        #inserts the correct values for the variable placeholders
        awg_program = awg_program.replace("FACT1", str(factor))
        awg_program = awg_program.replace("FACT2", str(factor2))
        awg_program = awg_program.replace("AMP1", str(amp))
        awg_program = awg_program.replace("AMP2", str(amp2))
        awg_program = awg_program.replace("P1", str(phase1))
        awg_program = awg_program.replace("P2", str(phase2))
        awg_program = awg_program.replace("PER", str(per))
        awg_program = awg_program.replace("WN1", str(WN1))
        awg_program = awg_program.replace("WN2", str(per2))
        awg_program = awg_program.replace("LEN1", str(round(LEN1,5))) #additional rounding is included to mitigate python's floating point errors.
        awg_program = awg_program.replace("LEN2", str(round(LEN2,5))) #additional rounding is included to mitigate python's floating point errors.
        
        #if the singletone is nested inside of a while loop
        if self.While==True: 
            self.seq= self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=awg_program

        if self.While==True:
            self.seq +="}"

        
    
    def addsingletoneRF(
        #This function writes a single-tone channel 1&2 signal+ a channel3 digital modulation 
        # based single-tone into the sequencer. 
        self,
        dur:float       = 0.001,
        freq:float      = 1e6,
        amp:float       = 1.0,
        amp2:float      = None,
        phase1:float    = None,
        phase2:float    = None,
        RFfreq:float    = None, #frequency of channel 3 output
        RFamp:float     = 1, #amplitude of channel 3 output
        RFphase:float   = 0, #absolute phase of channel 3 output (not mediated by phasetime)
        Mark:list       = [0,0,0,0],#if and where a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto phase1 and phase2.
    ):
        Flag1,Flag2=False,False
    
        if amp2==None:
            amp2=amp
            Flag1=True
        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
        self.t+=dur
        
        awg_program = textwrap.dedent(#initializes the string written to the sequencer
            """
            """
        )
        if RFfreq!=None: 
            #Checks if the frequency of channel 3 needs to be changed. If it does, a 20 ns break is inserted into the sequence
            # during which the program switches the frequency of the second oscillator, the phase of the 3rd channel,
            # and resets the oscillator phase to 0.
            dur-=2e-7
            awg_program+=textwrap.dedent("""
                playZero(240);
                setInt("oscs/1/freq", RFfreq);
                setDouble("sines/0/phaseshift",P1);
                playZero(240);
                resetOscPhase();""".replace("RFfreq",str(RFfreq)).replace("P1",str(RFphase*180/np.pi))
            )
            phasetime+=2e-7

        

        #default phases of channel 1 and 2 are sine and cosine.        
        if phase1==None:
            phase1=0
        if phase2==None:
            phase2=phase1-np.pi/2
            Flag2=True
        
        damp,dphase=HDAWGcalibration(amp,amp2,freq,phase1,phase2)
        if Flag1:
            amp2+=damp
        if Flag2:
            phase2+=dphase

        # print('bean')
        #For various calcs, need the following factors.
        per=freq*dur
        if dur>1e-7:
            # print(dur)
            factor=max(int(np.log2(2.4e9/(self.acc*freq))),0)
        else:
            factor=0
        
        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        if freq in self.freqdict:
            [WN1,LEN1,freqout,factor2]=self.freqdict[freq]
            if dur>1e-7 and factor!=factor2:
                WN1,LEN1,freqout=findfreq(freq,self.Ecc,self.msamp,factor)
                factor=factor2
        else:
            WN1,LEN1,freqout=findfreq(freq,self.Ecc,self.msamp,factor)
            self.freqdict[freq]=[WN1,LEN1,freqout,factor]

        N=int(round(2.4e9*dur*1/(2**factor),5))

        #Interpreting part of the output of freqout()
        if WN1==0:
            per=0
        else:
            per=int(per/WN1)
        
        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addsingletone does not)
        phase1+=(freqout*2*np.pi*phasetime)%(2*np.pi)
        phase2+=(freqout*2*np.pi*phasetime)%(2*np.pi)


        #Find the length of the final, unrepeated, waveform. Checking to ensure the length of this waveform is positive.
        LEN2=N-LEN1*per
        while(LEN2<0):
            per=per-1
            LEN2=N-LEN1*per
        per2=freqout*(LEN2/24e8)*2**factor

        #Makes sure the final waveform is long enough
        if round(LEN2,5)<100 and per>0:
            per-=1
            LEN2+=LEN1
            per2+=WN1
        
        factor2=factor
        while ((LEN2%16!=0) and (factor2>0)):
            factor2=factor2-1
            LEN2=2*LEN2
        if LEN2%16!=0:
            print("The precision required by this pulse is not achievable. Will zero extend this part of the pulse, which may cause various issues.")
        
        
        #Starting the awg_program for this section depending on whether there should be a repeating section:
        if per==0:
            awg_program += textwrap.dedent(
                """
                    playWave(1,marker(LEN2,1)+OUTP1,2,marker(LEN2,1)+OUTP2,3,marker(LEN2,1)+RFOUTP,4,marker(LEN2,1),FACT2);
                """
            )
        else:
            awg_program += textwrap.dedent(
                """
                    repeat(PER){
                        playWave(1,marker(LEN1,1)+OUT1,2,marker(LEN1,1)+OUT2,3,marker(LEN1,1)+RFOUTI,4,marker(LEN1,1),FACT1);
                    }
                    playWave(1,marker(LEN2,1)+OUTP1,2,marker(LEN2,1)+OUTP2,3,marker(LEN2,1)+RFOUTP,4,marker(LEN2,1),FACT2);
                """
            )
        
        
        #Ensures the markers are inserted into the correct channels, and placeholders are erased accordingly
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1)+",str(i+1)+",")
        
        #The following code replaces the placeholders with the correct waveforms
        OUT1="AMP1*sine(LEN1, P1, WN1)"
        OUT2="AMP2*sine(LEN1, P2, WN1)"
        RFOUTI="rect(LEN1,RFAMP)"
        awg_program=awg_program.replace("OUT1",OUT1)
        awg_program=awg_program.replace("OUT2",OUT2)
        awg_program=awg_program.replace("RFOUTI",RFOUTI)
        OUTP1="AMP1*sine(LEN2, P1,WN2)"
        OUTP2="AMP2*sine(LEN2, P2,WN2)"
        RFOUTP="rect(LEN2,RFAMP)"
        awg_program=awg_program.replace("OUTP1",OUTP1)
        awg_program=awg_program.replace("OUTP2",OUTP2)
        awg_program=awg_program.replace("RFOUTP",RFOUTP)

        #inserts the correct values for the variable placeholders
        awg_program=awg_program.replace("RFAMP",str(RFamp))
        awg_program = awg_program.replace("FACT1", str(factor))
        awg_program = awg_program.replace("FACT2", str(factor2))
        awg_program = awg_program.replace("AMP1", str(amp))
        awg_program = awg_program.replace("AMP2", str(amp2))
        awg_program = awg_program.replace("P1", str(phase1))
        awg_program = awg_program.replace("P2", str(phase2))
        awg_program = awg_program.replace("PER", str(per))
        awg_program = awg_program.replace("WN1", str(WN1))
        awg_program = awg_program.replace("WN2", str(per2))
        awg_program = awg_program.replace("LEN1", str(round(LEN1,5)))
        awg_program = awg_program.replace("LEN2", str(round(LEN2,5)))

        
        if self.While==True:#if the singletone is nested inside of a while loop
            self.seq= self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=awg_program
        if self.While==True:
            self.seq +="}"
        

    def addtwotoneRF(
        #Produces a two-tone two-channel signal using the awg-sequencer.
        self,#sequence to append
        dur:float = 0.0001, #duration
        freq1:float = 1e7, #frequency 1st tone
        amp11:float = 0.5, #amplitude 1st tone
        amp12:float = None, #amplitude 1st tone
        freq2:float = 3e6, #frequency 2nd tone
        amp21:float = 0.5, #amplitude 2nd tone
        amp22:float = None, #amplitude 2nd tone
        phase11:float=None, #phase 1st tone 1st channel
        phase12:float=None, #phase 1st tone 2st channel
        phase21:float=None, #phase 2st tone 1st channel
        phase22:float=None, #phase 2st tone 2st channel
        RFfreq:float    = None, #frequency of channel 3 output
        RFamp:float     = 1, #amplitude of channel 3 output
        RFphase:float   = 0, #absolute phase of channel 3 output (not mediated by phasetime)
        Mark:list=[0,0,0,0], #if a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto the phases.
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2  #2nd Output channel
    ):
        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
        self.t+=dur

        Flag11,Flag12,Flag21,Flag22=False,False,False,False
        #default phases of channel 1 and 2 are sine and cosine.        
        if phase11==None:
            phase11=0
        if phase12==None:
            phase12=phase11-np.pi/2
            Flag11=True
        
        #default phases of channel 1 and 2 are sine and cosine.        
        if phase21==None:
            phase21=0
        if phase22==None:
            phase22=phase21-np.pi/2
            Flag12=True
        
        #add method for channel dependant amplitudes
        if amp12==None:
            amp12=amp11
            Flag21=True
        #add method for channel dependant amplitudes
        if amp22==None:
            amp22=amp21
            Flag22=True

        awg_program = textwrap.dedent(#initializes the string written to the sequencer
            """
            """
        )
        if RFfreq!=None: 
            #Checks if the frequency of channel 3 needs to be changed. If it does, a 20 ns break is inserted into the sequence
            # during which the program switches the frequency of the second oscillator, the phase of the 3rd channel,
            # and resets the oscillator phase to 0.
            dur-=2e-7
            awg_program+=textwrap.dedent("""
                playZero(240);
                setInt("oscs/1/freq", RFfreq);
                setDouble("sines/0/phaseshift",P1);
                playZero(240);
                resetOscPhase();""".replace("RFfreq",str(RFfreq)).replace("P1",str(RFphase*180/np.pi))
            )
            phasetime+=2e-7

        #factor by which the sampling rate must be reduced (sampling rate=2.4e9/(2**factor))
        # factor=min(max(int(np.log2(2.4e9/(self.acc*max(freq1,freq2)))),0),1)
        
        if dur>1e-7:
            factor=min(max(int(np.log2(2.4e9/(self.acc*max(freq1,freq2)))),0),1)
        else:
            factor=0

        N=int(round(2.4e9*dur*1/(2**factor),5))
        while N%16!=0 and factor>0:
            factor-=1
            N=int(round(2.4e9*dur*1/(2**factor),5))
        if N%16!=0:
            print("The precision required by this pulse is not achievable. Will zero extend this part of the pulse, which will cause various issues.")
            

        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        #In the case of the two-tone no approximations are actually needed, but for the sake of consistency, the approximation is still used.
        if freq1 in self.freqdict:
            [_,_,freq1out,_]=self.freqdict[freq1]
        else:
            WN1,LEN1,freq1out=findfreq(freq1,self.Ecc,self.msamp,factor)
            self.freqdict[freq1]=[WN1,LEN1,freq1out,factor]

        if freq2 in self.freqdict:
            [_,_,freq2out,_]=self.freqdict[freq2]
        else:
            WN1,LEN1,freq2out=findfreq(freq2,self.Ecc,self.msamp,factor)
            self.freqdict[freq2]=[WN1,LEN1,freq2out,factor]

    
        #For various calcs, need the following constants.
        per1=freq1out*dur
        per2=freq2out*dur
        
                

        
        damp,dphase=HDAWGcalibration(amp11,amp12,freq1,phase11,phase12)
        if Flag21:
            amp12+=damp
        if Flag11:
            phase12+=dphase
        damp,dphase=HDAWGcalibration(amp21,amp22,freq2,phase21,phase22)
        if Flag22:
            amp22+=damp
        if Flag12:
            phase22+=dphase

        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addtwotone does not)
        phase11+=(freq1out*phasetime*np.pi*2)%(2*np.pi)
        phase12+=(freq1out*phasetime*np.pi*2)%(2*np.pi)
        phase21+=(freq2out*phasetime*np.pi*2)%(2*np.pi)
        phase22+=(freq2out*phasetime*np.pi*2)%(2*np.pi)
        
        #initializing code inserted into the sequencer, with all placeholders in place
        awg_program+=textwrap.dedent("""
        playWave(1,marker(Num,1)+WAVE1,2,marker(Num,1)+WAVE2,3,marker(Num,1)+WAVE3+RFOUT,4,marker(Num,1)+WAVE4,Fact);""")

        #Replacing the correct placeholders to output in the correct channels.
        OUT1="+add(Amp11*sine(Num,Phase11,Per1),Amp21*sine(Num,Phase21,Per2))"
        OUT2="+add(Amp12*sine(Num,Phase12,Per1),Amp22*sine(Num,Phase22,Per2))"
        if CHAN1==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT1)
        else:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT2)
        awg_program=awg_program.replace("+WAVE1","")
        awg_program=awg_program.replace("+WAVE2","")
        awg_program=awg_program.replace("+WAVE3","")
        awg_program=awg_program.replace("+WAVE4","")

        
        RFOUT="rect(Num,RFAMP)"
        awg_program=awg_program.replace("RFOUT",RFOUT)

        #Ensures the markers are inserted into the correct channels, and placeholders are erased accordingly
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(Num,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(Num,1),","")

        #replaces the placeholders with the correct values.
        awg_program=awg_program.replace("Num", str(round(N,5)))
        awg_program=awg_program.replace("Amp11", str(amp11))
        awg_program=awg_program.replace("Amp12", str(amp12))
        awg_program=awg_program.replace("Amp21", str(amp21))
        awg_program=awg_program.replace("Amp22", str(amp22))
        awg_program=awg_program.replace("Per1", str(per1))
        awg_program=awg_program.replace("Per2", str(per2))
        awg_program=awg_program.replace("Fact", str(factor))
        awg_program=awg_program.replace("Phase11", str(phase11))
        awg_program=awg_program.replace("Phase12", str(phase12))
        awg_program=awg_program.replace("Phase21", str(phase21))
        awg_program=awg_program.replace("Phase22", str(phase22))
        awg_program=awg_program.replace("RFAMP",str(RFamp))


        #if the singletone is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequuence
        self.seq+=awg_program
        if self.While==True:
            self.seq+="\n}"
    
    def modulatorsingletone(self,tdur=1e-3,freq=None,phase1=None,phase2=None,amp=1,amp2=None,Mark=[0,0,0,0]):
        #Writes code into the sequencer for a digital modulator based two-channel single-tone.
        #The code assumes that the digital modulators have been enabled.
        #No method for turning these on/off within the sequencer has been found yet.

        Flag1,Flag2=False,False

        #add method for channel dependant amplitudes
        if amp2==None:
            amp2=amp
            Flag1=True

        k=False
        if phase1==None:
            phase1=0
            k=True
        if phase2==None:
            phase2=phase1-np.pi/2
            Flag2=True
            k=True
        
        damp,dphase=HDAWGcalibration(amp,amp2,freq,phase1,phase2) #implements calibration
        
        if Flag1:
            amp2+=damp
        if Flag2:
            phase2+=dphase

        #initial values helpful with further calculations
        per=int(tdur*24e8/240)-1
        LEN1=240
        LEN2=int(tdur*24e8-per*240)
        
        #initializes the program
        awg_program = textwrap.dedent("""
        """)

        #If the frequency is changed, change the frequency.
        if freq!=None:
            if per>0:
                per-=1
            else:
                LEN2-=240
            awg_program+=textwrap.dedent("""
            playZero(176);
            setDouble("oscs/1/freq", %s);
            playZero(64);
            resetOscPhase();""" %str(freq)

            )
        
        #If the phase is changed, change it accordingly.
        if k:
            awg_program+=textwrap.dedent("""
            setInt("sines/2/phaseshift",P1);
            setInt("sines/3/phaseshift",P2);
            resetOscPhase();""".replace("P1",str(phase1*180/np.pi)).replace("P2",str(phase2*180/np.pi))
            )

        if per>0:#adds the code to produce the rectangular pulse needed for digital modulation, with a check to reduce redundant code
            awg_program += textwrap.dedent(
                """
                    repeat(PER){
                        playWave(1,marker(LEN1,1)+AMP1*ones(LEN1),2,marker(LEN1,1)+AMP2*ones(LEN1),3,marker(LEN1,1),4,marker(LEN1,1),0);
                    }
                    playWave(1,marker(LEN2,1)+AMP1*ones(LEN2),2,marker(LEN2,1)+AMP2*ones(LEN2),3,marker(LEN2,1),4,marker(LEN2,1));
                """
            )
        else:
            awg_program += textwrap.dedent(
                """
                    playWave(1,marker(LEN2,1)+AMP1*ones(LEN2),2,marker(LEN2,1)+AMP2*ones(LEN2),3,marker(LEN2,1),4,marker(LEN2,1));

                """
            )

        #Adds the markers, and removes placeholders
        for i in range(4):
            if Mark[i]==0:
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN1,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1),","")
                awg_program=awg_program.replace(str(i+1)+",marker(LEN2,1)+",str(i+1)+",")
                
        #inserts the correct values for the variable placeholders
        awg_program = awg_program.replace("PER", str(per))
        awg_program = awg_program.replace("LEN1", str(round(LEN1,5)))
        awg_program = awg_program.replace("LEN2", str(round(LEN2,5)))
        awg_program = awg_program.replace("AMP1", str(amp))
        awg_program = awg_program.replace("AMP2", str(amp2))

        #if the singletone is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequence
        self.seq+=awg_program
        if self.While==True:
            self.seq+="\n}"

    def addtwotone(
        #Produces a two-tone two-channel signal using the awg-sequencer.
        self,#sequence to append
        dur:float = 0.0001, #duration
        freq1:float = 1e7, #frequency 1st tone
        amp11:float = 0.5, #amplitude 1st tone
        amp12:float = None, #amplitude 1st tone
        freq2:float = 3e6, #frequency 2nd tone
        amp21:float = 0.5, #amplitude 2nd tone
        amp22:float = None, #amplitude 2nd tone
        phase11:float=None, #phase 1st tone 1st channel
        phase12:float=None, #phase 1st tone 2st channel
        phase21:float=None, #phase 2st tone 1st channel
        phase22:float=None, #phase 2st tone 2st channel
        Mark:list=[0,0,0,0], #if a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto the phases.
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2  #2nd Output channel
    ):
        #automatically assigns phasetime if it is left undefined
        if phasetime==None:
            phasetime=self.t
        self.t+=dur

        Flag11,Flag12,Flag21,Flag22=False,False,False,False
        #default phases of channel 1 and 2 are sine and cosine.        
        if phase11==None:
            phase11=0
        if phase12==None:
            phase12=phase11-np.pi/2
            Flag11=True
        
        #default phases of channel 1 and 2 are sine and cosine.        
        if phase21==None:
            phase21=0
        if phase22==None:
            phase22=phase21-np.pi/2
            Flag12=True
        
        #add method for channel dependant amplitudes
        if amp12==None:
            amp12=amp11
            Flag21=True
        #add method for channel dependant amplitudes
        if amp22==None:
            amp22=amp21
            Flag22=True

        #factor by which the sampling rate must be reduced (sampling rate=2.4e9/(2**factor))
        
        if dur>1e-7:
            factor=min(max(int(np.log2(2.4e9/(self.acc*max(freq1,freq2)))),0),1)
        else:
            factor=0

        N=int(round(2.4e9*dur*1/(2**factor),5))
        while N%16!=0 and factor>0:
            factor-=1
            N=int(round(2.4e9*dur*1/(2**factor),5))
        if N%16!=0:
            print("The precision required by this pulse is not achievable. Will zero extend this part of the pulse, which will cause various issues.")
            
        #If this frequency has already been calculated for, use the previous results, else store the calculation result in the freqdict.
        #In the case of the two-tone no approximations are actually needed, but for the sake of consistency, the approximation is still used. The approximation should be good enough.
        if freq1 in self.freqdict:
            [_,_,freq1out,_]=self.freqdict[freq1]
        else:
            WN1,LEN1,freq1out=findfreq(freq1,self.Ecc,self.msamp,factor)
            self.freqdict[freq1]=[WN1,LEN1,freq1out,factor]
        if freq2 in self.freqdict:
            [_,_,freq2out,_]=self.freqdict[freq2]
        else:
            WN1,LEN1,freq2out=findfreq(freq2,self.Ecc,self.msamp,factor)
            self.freqdict[freq2]=[WN1,LEN1,freq2out,factor]

    
        #For various calcs, need the following constants.
        per1=freq1out*dur
        per2=freq2out*dur
        

        #Implements the Calibration.
        damp,dphase=HDAWGcalibration(amp11,amp12,freq1,phase11,phase12)
        if Flag21:
            amp12+=damp
        if Flag11:
            phase12+=dphase
        damp,dphase=HDAWGcalibration(amp21,amp22,freq2,phase21,phase22)
        if Flag22:
            amp22+=damp
        if Flag12:
            phase22+=dphase

        #define phase in terms of a time t0. (This phasetime needs to be assigned. Sequence.sequence_generator takes care of this automatically, but Sequence.addtwotone does not)
        phase11+=(freq1out*phasetime*np.pi*2)%(2*np.pi)
        phase12+=(freq1out*phasetime*np.pi*2)%(2*np.pi)
        phase21+=(freq2out*phasetime*np.pi*2)%(2*np.pi)
        phase22+=(freq2out*phasetime*np.pi*2)%(2*np.pi)
        
        #initializing code inserted into the sequencer, with all placeholders in place
        awg_program=textwrap.dedent("""
        playWave(1,marker(Num,1)+WAVE1,2,marker(Num,1)+WAVE2,3,marker(Num,1)+WAVE3,4,marker(Num,1)+WAVE4,Fact);""")

        #Replacing the correct placeholders to output in the correct channels.
        OUT1="+add(Amp11*sine(Num,Phase11,Per1),Amp21*sine(Num,Phase21,Per2))"
        OUT2="+add(Amp12*sine(Num,Phase12,Per1),Amp22*sine(Num,Phase22,Per2))"
        if CHAN1==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT1)
        elif CHAN1==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT1)
        elif CHAN1==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT1)
        else:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT1)
        if CHAN2==1:
            awg_program=awg_program.replace("WAVE1","WAVE1"+OUT2)
        if CHAN2==2:
            awg_program=awg_program.replace("WAVE2","WAVE2"+OUT2)
        if CHAN2==3:
            awg_program=awg_program.replace("WAVE3","WAVE3"+OUT2)
        if CHAN2==4:
            awg_program=awg_program.replace("WAVE4","WAVE4"+OUT2)
        awg_program=awg_program.replace("+WAVE1","")
        awg_program=awg_program.replace("+WAVE2","")
        awg_program=awg_program.replace("+WAVE3","")
        awg_program=awg_program.replace("+WAVE4","")

        #Ensures the markers are inserted into the correct channels, and placeholders are erased accordingly
        for i in range(4):
            if round(Mark[i])==0:
                awg_program=awg_program.replace(str(i+1)+",marker(Num,1)+",str(i+1)+",")
                awg_program=awg_program.replace(str(i+1)+",marker(Num,1),","")

        #replaces the placeholders with the correct values.
        awg_program=awg_program.replace("Num", str(round(N,5)))
        awg_program=awg_program.replace("Amp11", str(amp11))
        awg_program=awg_program.replace("Amp12", str(amp12))
        awg_program=awg_program.replace("Amp21", str(amp21))
        awg_program=awg_program.replace("Amp22", str(amp22))
        awg_program=awg_program.replace("Per1", str(per1))
        awg_program=awg_program.replace("Per2", str(per2))
        awg_program=awg_program.replace("Fact", str(factor))
        awg_program=awg_program.replace("Phase11", str(phase11))
        awg_program=awg_program.replace("Phase12", str(phase12))
        awg_program=awg_program.replace("Phase21", str(phase21))
        awg_program=awg_program.replace("Phase22", str(phase22))


        #if the singletone is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        #adding the code, now written, into the sequuence
        self.seq+=awg_program
        if self.While==True:
            self.seq+="\n}"

    def addtwotones(
        #Produces a two-tone two-channel signal using the awg-sequencer. This particular bit of code splits up a two-tone signal into multiple two tones in order to conserve cache memory.
        #The default timesplit has been taken to be 75 us as this is approximately 50% of the cache, allowing two splits to be loaded into the cache simultaneously.
        self,#sequence to append
        dur:float = 0.0001, #duration
        freq1:float = 1e7, #frequency 1st tone
        amp11:float = 0.5, #amplitude 1st tone
        amp12:float = None, #amplitude 1st tone
        freq2:float = 3e6, #frequency 2nd tone
        amp21:float = 0.5, #amplitude 2nd tone
        amp22:float = None, #amplitude 2nd tone
        phase11:float=None, #phase 1st tone 1st channel
        phase12:float=None, #phase 1st tone 2st channel
        phase21:float=None, #phase 2st tone 1st channel
        phase22:float=None, #phase 2st tone 2st channel
        Mark:list=[0,0,0,0], #if a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto the phases.
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2,  #2nd Output channel
        timesplit:float = 80e-6,
    ):
        dur=round(dur*2.4e9,4)/2.4e9
        if phasetime==None:
            phasetime=0
        if dur>timesplit:
            for i in range(int(dur/timesplit)):
                t=phasetime+i*timesplit
                self.addtwotone(dur=timesplit,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,Mark=Mark,phasetime=t,CHAN1=CHAN1,CHAN2=CHAN2)
            if round(2.4e9*(dur/timesplit),5)>timesplit*2.4e9*(i+1)+240:
                t+=timesplit
                self.addtwotone(dur=dur%timesplit,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,Mark=Mark,phasetime=t,CHAN1=CHAN1,CHAN2=CHAN2)
        else:
            self.addtwotone(dur=dur,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,Mark=Mark,phasetime=phasetime,CHAN1=CHAN1,CHAN2=CHAN2)


    def addtwotonesRF(
        #Produces a two-tone two-channel signal using the awg-sequencer.
        self,#sequence to append
        dur:float = 0.0001, #duration
        freq1:float = 1e7, #frequency 1st tone
        amp11:float = 0.5, #amplitude 1st tone
        amp12:float = None, #amplitude 1st tone
        freq2:float = 3e6, #frequency 2nd tone
        amp21:float = 0.5, #amplitude 2nd tone
        amp22:float = None, #amplitude 2nd tone
        phase11:float=None, #phase 1st tone 1st channel
        phase12:float=None, #phase 1st tone 2st channel
        phase21:float=None, #phase 2st tone 1st channel
        phase22:float=None, #phase 2st tone 2st channel
        Mark:list=[0,0,0,0], #if a marker is desired
        phasetime:float = None, #Time evolution of the phase, this time*freq_real is added onto the phases.
        RFfreq:float    = None, #frequency of channel 3 output
        RFamp:float     = 1, #amplitude of channel 3 output
        RFphase:float   = 0, #absolute phase of channel 3 output (not mediated by phasetime)
        CHAN1:int       = 1, #1st Output channel
        CHAN2:int       = 2,  #2nd Output channel
        timesplit:float = 80e-6,
    ):
        dur=round(dur*2.4e9,4)/2.4e9
        if phasetime==None:
            phasetime=0
        if dur>(timesplit+1):
            self.addtwotoneRF(dur=timesplit,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,RFfreq=RFfreq,RFamp=RFamp,RFphase=RFphase,Mark=Mark,phasetime=phasetime,CHAN1=CHAN1,CHAN2=CHAN2)
            phasetime+=timesplit
            for i in range(int(dur/timesplit)-1):
                t=phasetime+i*timesplit
                self.addtwotoneRF(dur=timesplit,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,RFfreq=RFfreq,RFamp=RFamp,RFphase=RFphase,Mark=Mark,phasetime=phasetime,CHAN1=CHAN1,CHAN2=CHAN2)
            if round(2.4e9*(dur/timesplit),5)>timesplit*2.4e9*(i+1)+240:
                t+=timesplit
                self.addtwotoneRF(dur=dur%timesplit,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,RFfreq=RFfreq,RFamp=RFamp,RFphase=RFphase,Mark=Mark,phasetime=phasetime,CHAN1=CHAN1,CHAN2=CHAN2)
        else:
            self.addtwotoneRF(dur=dur,freq1=freq1,amp11=amp11,amp12=amp12,freq2=freq2,amp21=amp21,amp22=amp22,phase11=phase11,phase12=phase12,phase21=phase21,phase22=phase22,RFfreq=RFfreq,RFamp=RFamp,RFphase=RFphase,Mark=Mark,phasetime=phasetime,CHAN1=CHAN1,CHAN2=CHAN2)
            

    def addwait(
        #adds a wait to the sequence, with or without markers
        self,
        dur:float = 0.001,#duration of wait
        Mark:list=[0,0,0,0],#if and where a marker is desired
        
    ):
        #adds onto time for phasetime calculations.
        self.t+=dur
        #if the wait is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        

        if Mark==[0,0,0,0]: #if no markers are wanted then, the playzero command is used as it takes up close to 0 memory
            self.seq+=textwrap.dedent("""
                repeat(B){
                    playZero(240000000);
                }
            playZero(D);""").replace("D",str(round((dur%0.1)*24e8))).replace("B",str(int(dur*10)))

            
        else:# if there are markers, then there must be a playWave function, which can then be made up of repeated short fragments.
            self.seq+=textwrap.dedent("""
            repeat(NREP){
                playWave(1,marker(240,AUX1),2,marker(240,AUX2),3,marker(240,AUX3),4,marker(240,AUX4));
            }
            playWave(1,marker(N2nd,AUX1),2,marker(N2nd,AUX2),3,marker(N2nd,AUX3),4,marker(N2nd,AUX4));""")

            NREP=int(round(dur*1e7,8))-1
            if NREP==-1:
                NREP=0
                N2nd=int(round(dur*24e8,5))
            else:
                N2nd=int(round(dur*24e8,5))%240+240

            self.seq=self.seq.replace("NREP",str(NREP))
            self.seq=self.seq.replace("N2nd",str(N2nd))
            #A simpler, more memory-intensive method of assigning markers is used here. This is unlikely to cause issues due to the shortness of these fragments.
            self.seq=self.seq.replace("AUX1",str(Mark[0]))
            self.seq=self.seq.replace("AUX2",str(Mark[1]))
            self.seq=self.seq.replace("AUX3",str(Mark[2]))
            self.seq=self.seq.replace("AUX4",str(Mark[3]))
        if self.While==True:
            self.seq+="\n}"
        
        
    def phasch(self,
               RFphase=0):
        self.seq+=textwrap.dedent("""
                setDouble("sines/0/phaseshift",P1);
                """.replace("P1",str(RFphase*180/np.pi))
            )
        
        
    def addwaitRF(
        #adds a wait to the sequence, with or without markers
        self,
        dur:float = 0.001,#duration of wait
        RFfreq:float    = None, #frequency of channel 3 output
        RFamp:float     = 1, #amplitude of channel 3 output
        RFphase:float   = 0, #absolute phase of channel 3 output (not mediated by phasetime)
        Mark:list=[0,0,0,0],#if and where a marker is desired
        
    ):
        #adds onto time for phasetime calculations.
        self.t+=dur
        #if the wait is nested inside of a while loop
        if self.While==True:
            self.seq=self.seq[:len(self.seq)-2]
        
        if RFfreq!=None: 
            #Checks if the frequency of channel 3 needs to be changed. If it does, a 20 ns break is inserted into the sequence
            # during which the program switches the frequency of the second oscillator, the phase of the 3rd channel,
            # and resets the oscillator phase to 0.
            dur-=2e-7
            self.seq+=textwrap.dedent("""
                playZero(240);
                setInt("oscs/1/freq", RFfreq);
                setDouble("sines/0/phaseshift",P1);
                playZero(240);
                resetOscPhase();""".replace("RFfreq",str(RFfreq)).replace("P1",str(RFphase*180/np.pi))
            )

        
        self.seq+=textwrap.dedent("""
        repeat(NREP){
            playWave(1,marker(240,AUX1),2,marker(240,AUX2),3,marker(240,AUX3)+rect(240,RFAMP),4,marker(240,AUX4));
        }
        playWave(1,marker(N2nd,AUX1),2,marker(N2nd,AUX2),3,marker(N2nd,AUX3),4,marker(N2nd,AUX4));""")

        NREP=int(round(dur*1e7,8))-1
        if NREP==-1:
            NREP=0
            N2nd=int(round(dur*24e8,5))
        else:
            N2nd=int(round(dur*24e8,5))%240+240

        self.seq=self.seq.replace("NREP",str(NREP))
        self.seq=self.seq.replace("N2nd",str(N2nd))
        #A simpler, more memory-intensive method of assigning markers is used here. This is unlikely to cause issues due to the shortness of these fragments.
        self.seq=self.seq.replace("AUX1",str(Mark[0]))
        self.seq=self.seq.replace("AUX2",str(Mark[1]))
        self.seq=self.seq.replace("AUX3",str(Mark[2]))
        self.seq=self.seq.replace("AUX4",str(Mark[3]))
        self.seq=self.seq.replace("RFAMP",str(RFamp))
        if self.While==True:
            self.seq+="\n}"
    



    
            


    def Upload(self,device,daq):
        #uploads the sequence to the HDAWG.
        #This code is ?s?t?o?l?e?n?  borrowed from the Zurich Instruments example code. I am not sure how it works, but it has not failed me yet.


        # Create an instance of the AWG Module
        awgModule = daq.awgModule()
        awgModule.set("device", device)
        awgModule.execute()

        # Get the modules data directory
        data_dir = awgModule.getString("directory")
        # All CSV files within the waves directory are automatically recognized by the AWG module
        wave_dir = path.join(data_dir, "awg", "waves")
        if not path.isdir(wave_dir):
            # The data directory is created by the AWG module and should always exist. If this exception
            # is raised, something might be wrong with the file system.
            raise Exception(
                f"AWG module wave directory {wave_dir} does not exist or is not a directory"
            )

        # Save waveform data to CSV
        # Transfer the AWG sequence program. Compilation starts automatically.
        awgModule.set("compiler/sourcestring", self.seq)
        # Note: when using an AWG program from a source file (and only then), the compiler needs to
        # be started explicitly with awgModule.set('compiler/start', 1)
        while awgModule.getInt("compiler/status") == -1:
            sleep(0.1)

        if awgModule.getInt("compiler/status") == 1:
            # compilation failed, raise an exception
            raise Exception(awgModule.getString("compiler/statusstring"))

        if awgModule.getInt("compiler/status") == 0:
            print(
                "Compilation successful with no warnings, will upload the program to the instrument."
            )
        if awgModule.getInt("compiler/status") == 2:
            print(
                "Compilation successful with warnings, will upload the program to the instrument."
            )
            print("Compiler warning: ", awgModule.getString("compiler/statusstring"))
            # print(self.seq)
            raise RuntimeError("HDAWG compilation warning, see lib/site-pacakges/labscript_devices/HDAWG/HDAWG_funcs line 1233")
        i = 0
        while (awgModule.getDouble("progress") < 1.0) and (
            awgModule.getInt("elf/status") != 1
        ):
            print(f"{i} progress: {awgModule.getDouble('progress'):.2f}", end="\r")
            #time.sleep(0.2)
            i += 1
        print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
        if awgModule.getInt("elf/status") == 0:
            print("Upload to the instrument successful.")
        if awgModule.getInt("elf/status") == 1:
            raise Exception("Upload to the instrument failed.")
    


        # This is the preferred method of using the AWG: Run in single mode continuous waveform playback
        # is best achieved by using an infinite loop (e.g., while (true)) in the sequencer program.
        daq.setInt(f"/{device}/awgs/0/single", 1)
        daq.setInt(f"/{device}/awgs/0/enable", 1)


        print("Sequence Activated, Waiting on trigger")
    


    
    def sequence_generator(self,Pulses,Markers):
        #This function takes in a set of pulses and markers given by an external source and turns it into a sequence.
        #This code is deprecated and incompatible with the methodology used in labscript.
        Pulses.sort(key=operator.itemgetter("t0"))


        timesp,timesm,tymesp,tymesm=[],[],[],[]
        for P in Pulses:
            #Creates arrays of timestamps for the pulses
            timesp.append(P["t0"])
            if P["t0"]<0:
                raise(Exception("Pulses must occur after the trigger"))
            timesm.append(P["t0"]+P["tdur"])
            if P["tdur"]<=0:
                raise(Exception("Pulses must have positive duration"))
        
        
        for M in Markers:
            #Creates arrays of timestamps for the markers
            tymesp.append(M["t0"])
            if M["t0"]<0:
                raise(Exception("Markers must occur after the trigger"))
            tymesm.append(M["t0"]+M["tdur"])
            if M["tdur"]<=0:
                raise(Exception("Pulses must have positive duration"))

        #Checks whether valid timestamps were given
        for m in range(len(tymesp)):
            if m<len(tymesp)-1:
                if tymesp[m+1]-tymesm[m]<0:
                    raise(Exception("Marker",str(m+1),"and",str(m+2),"overlap. This is not allowed"))

        for t in range(len(timesp)):
            if t<len(timesp)-1:
                if timesp[t+1]-timesm[t]<0:
                    raise(Exception("Pulses",str(t+1),"and",str(t+2),"overlap. This is not allowed"))
        
        #Creates a full list for any changes to the sequence
        Events=[*set(tymesp+tymesm+timesm+timesp)]
        Events.sort()
        Mark=[0,0,0,0]
        
        #if the sequence should wait after triggering
        if Events[0]>0:
            self.addwait(dur=Events[0],Mark=False)
        P=0
        #Based on the timestamps and events, generate the pulse sequence.
        for m in range(len(Events)-1):
            t=Events[m]
            dur=Events[m+1]-t
            if t in tymesp:
                n=tymesp.index(t)
                Mark=Markers[n]["AUX"]
            elif t in tymesm:
                Mark=[0,0,0,0]

            if t in timesp:
                n=timesp.index(t)
                P=Pulses[n]
            elif t in timesm:
                P=0
            
            

            if P==0:
                self.addwait(dur=dur,Mark=Mark)
            else:
                if P["type"]=="Single":
                    freq=P["freq"]
                    P1=P["phase1"]
                    P2=P["phase2"]
                    self.addsingletone(dur=dur,freq=freq,amp=P["Amp"],phase1=P1,phase2=P2,Mark=Mark,phasetime=t)
                elif P["type"]=="Double":
                    freq1=P["freq1"]
                    freq2=P["freq2"]
                    P11=P["phase11"]
                    P12=P["phase12"]
                    P21=P["phase21"]
                    P22=P["phase22"]
                    self.addtwotone(dur=dur,freq1=freq1,amp1=P["Amp1"],freq2=freq2,amp2=P["Amp2"],phase11=P11,phase12=P12,phase21=P21,phase22=P22,Mark=Mark)
            





def add_pulse(                          #Function for generating the pulse table used by the sequence_generator
    t0:float=0,                         #Start time in seconds
    tdur:float=2e-4,                    #Duration in seconds
    Type:str="Single",                  #Single/two-tone
    amp1:float=1.0,                     #Amplitude, or Amplitude first tone
    amp2:float=1.0,                     #Amplitude second tone
    freq1:float=1e7,                    #Frequency, or Frequency first tone
    freq2:float=2e7,                    #Frequency second tone
    phase11:float=0,                    #Phase offset first channel first tone
    phase12:float=np.pi/2,              #Phase offset second channel first tone
    phase21:float=0,                    #Phase offset first channel second tone
    phase22:float=np.pi/2,              #Phase offset second channel second tone

):
    P={}
    P["t0"]=t0
    P["tdur"]=tdur
    P["type"]=Type
    if Type=="Single": 
        P["Amp"]=amp1
        P["freq"]=freq1
        P["phase1"]=phase11
        P["phase2"]=phase12
    elif Type=="Double": 
        P["Amp1"]=amp1
        P["freq1"]=freq1
        P["phase11"]=phase11
        P["phase12"]=phase12
        P["Amp2"]=amp2
        P["freq2"]=freq2
        P["phase21"]=phase21
        P["phase22"]=phase22
    else:
        raise(Exception("Type of pulse not Recognized"))
    return P

def add_marker(                         #Function for generating the Marker table used by the sequence_generator
    t0:float=0,                         #Start time in seconds best results if 
    tdur:float=2e-4,                    #Duration in seconds
    AUXset:list=[1,0,0,0],
):
    M={}
    M["t0"]=t0
    M["tdur"]=tdur
    M["AUX"]=AUXset
    return M

    