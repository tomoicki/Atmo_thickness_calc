import csv
import io
import math
import numpy
import tkinter as tk
from ast import literal_eval
from tkinter import ttk

def is_num(x): #function which will check if input was a number
    try:
        float(x)
        return True
    except ValueError:
        return False

reader = csv.reader(open('Locations.txt', 'r'), delimiter=':')

Locations = []
for line in reader:
    Locations.append(line)
print(Locations[2],type(Locations[2]))
root = tk.Tk() 

#Date
DateFrame = tk.Frame(root,
                     bd=1,
                     padx=10,
                     pady=10,
                     relief=tk.SUNKEN)
DateFrame.pack()
DFLabel = tk.Label(DateFrame,
                   text="Date: ").pack(side=tk.LEFT)
YearEntry = tk.Entry(DateFrame,
                     width=6,
                     justify='center')
YearEntry.pack(side=tk.LEFT)
YearEntry.insert(0,"Year")
def YearEntryClear(event):
    YearEntry.delete(0,"end")
YearEntry.bind("<Button-1>", YearEntryClear)
def YearEntryFill(event):
    YearEntry.insert(0,"Year")
def YearEntryOut(x):
    i=YearEntry.get()
    if is_num(i) == False:
        YearEntryClear(x)
        YearEntryFill(x)
YearEntry.bind("<FocusOut>",YearEntryOut)
mon = tk.StringVar()
MonthList = ["January","February","March","April","May","June","July","August","September","October","November","December"]
Month31List = ["January","March","May","July","August","October","December"]
Month30List = ["April","June","September","November"]
ComboMonthsDefault = tk.StringVar(root)
ComboMonthsDefault.set("Month")
ComboMonths = ttk.Combobox(DateFrame,
                           width=10,
                           textvariable=ComboMonthsDefault,
                           values=MonthList,
                           state="readonly",
                           justify='center')
ComboMonths.pack(side=tk.LEFT)
DayList = [1]*31
Day30List = [1]*31
Day29List = [1]*31
for i in range(31):
    DayList[i]=i+1
    Day30List[i]=i+1
    Day29List[i]=i+1
Day30List.remove(31)
del Day29List[29:31]
ComboDaysDefault = tk.StringVar(root)
ComboDaysDefault.set("Day")
ComboDays = ttk.Combobox(DateFrame,
                         width=4,
                         textvariable=ComboDaysDefault,
                         state=tk.DISABLED,
                         justify='center')
ComboDays.pack(side=tk.RIGHT)
def DayUnlock(event):
    ComboDays.configure(state="readonly")
ComboMonths.bind("<Button-1>",DayUnlock)
def NumberOfDays(event):
    m = ComboMonths.get()
    if m in Month31List:
        ComboDays.configure(values=DayList)
    if m in Month30List:
        ComboDays.configure(values=Day30List)
    if m == 'February':
        ComboDays.configure(values=Day29List)
ComboDays.bind("<Enter>",NumberOfDays)
#Location
LocationFrame = tk.Frame(root,
                         bd=1,
                         padx=10,
                         pady=10,
                         relief=tk.SUNKEN)
LocationFrame.pack()
LocationLabel = tk.Label(LocationFrame,
                         text="Location").pack(side=tk.TOP)
LocationOwnFrame = tk.Frame(LocationFrame,
                            bd=3,
                            relief=tk.RIDGE)
LocationOwnFrame.pack(side=tk.LEFT)
RBValue = tk.StringVar()
RBOwn = tk.Radiobutton(LocationOwnFrame,
                       text="Your Own",
                       variable=RBValue,
                       tristatevalue=0,
                       value=1)
RBOwn.grid(row=0,column=1)
LatitudeLabel = tk.Label(LocationOwnFrame,
                         text="Latitude [\N{DEGREE SIGN}]").grid(row=1,column=0)
LatitudeEntry = tk.Entry(LocationOwnFrame,
                         width=6,
                         justify='center',
                         state=tk.DISABLED)
LatitudeEntry.grid(row=2,column=0)
def LatitudeEntryClear(event):
    LatitudeEntry.delete(0,"end")
LatitudeEntry.bind("<Button-1>", LatitudeEntryClear)
def LatitudeEntryFill(event):
    LatitudeEntry.insert(0,"0.00")
def LatitudeEntryOut(x):
    l=LatitudeEntry.get()
    if is_num(l) == False or l == '0.00' or int(l) > 90 or int(l) < -90:
        LatitudeEntryClear(x)
        LatitudeEntryFill(x)
LatitudeEntry.bind("<FocusOut>",LatitudeEntryOut)
LongitudeLabel = tk.Label(LocationOwnFrame,
                         text="Longitude [\N{DEGREE SIGN}]").grid(row=1,column=1)
LongitudeEntry = tk.Entry(LocationOwnFrame,
                          width=7,
                          justify='center')
LongitudeEntry.grid(row=2,column=1)
LongitudeEntry.configure(state=tk.DISABLED)
def LongitudeEntryClear(event):
    LongitudeEntry.delete(0,"end")
LongitudeEntry.bind("<Button-1>", LongitudeEntryClear)
def LongitudeEntryFill(event):
    LongitudeEntry.insert(0,"0.00")
def LongitudeEntryOut(x):
    l=LongitudeEntry.get()
    if is_num(l) == False or l == '0.00' or int(l) > 180 or int(l) < -180:
        LongitudeEntryClear(x)
        LongitudeEntryFill(x)
LongitudeEntry.bind("<FocusOut>",LongitudeEntryOut)
TimeZoneLabel = tk.Label(LocationOwnFrame,
                         text="TimeZone").grid(row=1,column=2)
TimeZoneEntry = tk.Entry(LocationOwnFrame,
                         justify='center',
                         width=4,
                         state=tk.DISABLED)
TimeZoneEntry.grid(row=2,column=2)
def TimeZoneEntryClear(event):
    TimeZoneEntry.delete(0,"end")
TimeZoneEntry.bind("<Button-1>", TimeZoneEntryClear)
def TimeZoneEntryFill(event):
    TimeZoneEntry.insert(0,"0")
def TimeZoneEntryOut(x):
    tz=TimeZoneEntry.get()
    if is_num(tz) == False:
        TimeZoneEntryClear(x)
        TimeZoneEntryFill(x)
TimeZoneEntry.bind("<FocusOut>",TimeZoneEntryOut)
LocationPredefFrame = tk.Frame(LocationFrame,
                               bd=3,
                               relief=tk.RIDGE)
LocationPredefFrame.pack(side=tk.LEFT)
RBDef = tk.Radiobutton(LocationPredefFrame,
                       text="Predefined",
                       variable=RBValue,
                       tristatevalue=0,
                       value=2)
RBDef.pack()

ComboLocationsDefault = tk.StringVar(root)
CLV = tk.StringVar()
ComboLocationsDefault.set("Location | +N/-S | +E/-W | \u00B1GMT")
ComboLocations = ttk.Combobox(LocationPredefFrame,
                              width=50,
                              textvariable=ComboLocationsDefault,
                              values=Locations,
                              justify='center',
                              state=tk.DISABLED)
ComboLocations.pack(side=tk.RIGHT)
def LocationOwnUnlock(event):
    ComboLocations.configure(state=tk.DISABLED)
    LatitudeEntry.configure(state=tk.NORMAL)
    LongitudeEntry.configure(state=tk.NORMAL)
    TimeZoneEntry.configure(state=tk.NORMAL)
    LatitudeEntry.delete(0,"end")
    LatitudeEntry.insert(0,"0.00")
    LongitudeEntry.delete(0,"end")
    LongitudeEntry.insert(0,"0.00")
    TimeZoneEntry.delete(0,"end")
    TimeZoneEntry.insert(0,"0")
RBOwn.bind("<Button-1>", LocationOwnUnlock)
def LocationPredefUnlock(event):
    LatitudeEntry.delete(0,"end")
    LongitudeEntry.delete(0,"end")
    TimeZoneEntry.delete(0,"end")
    LatitudeEntry.configure(state=tk.DISABLED)
    LongitudeEntry.configure(state=tk.DISABLED)
    TimeZoneEntry.configure(state=tk.DISABLED)
    ComboLocations.configure(state=tk.NORMAL)
    ComboLocations.configure(state="readonly")
RBDef.bind("<Button-1>", LocationPredefUnlock)

#Atmosphere Thickness
AtmoFrame = tk.Frame(root,
                     bd=1,
                     padx=10,
                     pady=10,
                     relief=tk.SUNKEN)
AtmoFrame.pack()
AtmoLabel = tk.Label(AtmoFrame,
                     text="Atmosphere Thickness").grid(row=0,column=1)
AtmoFramePredef = tk.Frame(AtmoFrame,
                           bd=3,
                           relief=tk.RIDGE)
AtmoFramePredef.grid(row=1,column=0)
var15 = tk.IntVar()
var35 = tk.IntVar()
var55 = tk.IntVar()
var100 = tk.IntVar()
tk.Checkbutton(AtmoFramePredef,
               text="15km",
               indicatoron=0,
               variable=var15).pack(side=tk.LEFT)
tk.Checkbutton(AtmoFramePredef,
               text="35km",
               indicatoron=0,
               variable=var35).pack(side=tk.LEFT)
tk.Checkbutton(AtmoFramePredef,
               text="55km",
               indicatoron=0,
               variable=var55).pack(side=tk.LEFT)
tk.Checkbutton(AtmoFramePredef,
               text="100km",
               indicatoron=0,
               variable=var100).pack(side=tk.LEFT)
AtmoFrameInfo = tk.Frame(AtmoFrame)
AtmoFrameInfo.grid(row=1,column=1)
AtmoInfo = """\n\nTroposphere spans from Earth's surface (with average density = 1.2985 kg/m3) 
up to 15km high (with average density 0.3639 kg/m3 at Tropopause).
Stratosphere spans from 15km up to 55km (with average density 0.0020 kg/m3 at Stratopause).
Ozone Layer is a part of Stratopshere, it is responsible for absorption of UV radiation.
It ends around 35th kilometer (with average density  = 0.0132 kg/m3)
At the altitude of 100km lays Karman Line which defines the border between Earth's atmosphere and outer space.
That's why we have predefinded atmoshpere thinknesses = 15, 35, 55 and 100 km.

As you can see, atmosphere density decreases very rapidly with altitude. 
Therefore, if you want to use this program for sunbathing (checking how "Sun strength" changes with time of the day), 
we advise to use "35km".\n"""
def info():
    DisplayText.insert(tk.END, AtmoInfo, 'justcenter')
    DisplayText.yview_pickplace("end")
InfoButton = tk.Button(AtmoFrameInfo,
                       text="Info",
                       bd=3,
                       command=info,
                       relief=tk.GROOVE).pack()
AtmoFrameOwn = tk.Frame(AtmoFrame,
                        bd=3,
                        relief=tk.RIDGE)
AtmoFrameOwn.grid(row=1,column=2)
AtmoOwnLabel = tk.Label(AtmoFrameOwn,
                        text="Your Own").pack(side=tk.LEFT)
def OwnAtmoEntryUnlock(event):
    if CheckVar.get() == 0:
        OwnAtmoEntry.configure(state=tk.NORMAL)
        OwnAtmoEntry.insert(0,"0")
    if CheckVar.get() == 1:
        OwnAtmoEntry.delete(0,"end")
        OwnAtmoEntry.configure(state=tk.DISABLED)
CheckVar = tk.IntVar()
AtmoCheck = tk.Checkbutton(AtmoFrameOwn,
                           variable=CheckVar,
                           onvalue=1,
                           offvalue=0)
AtmoCheck.pack(side=tk.LEFT)
AtmoCheck.bind("<Button-1>", OwnAtmoEntryUnlock)
OwnAtmoEntry = tk.Entry(AtmoFrameOwn,
                        width=6,
                        state=tk.DISABLED,
                        justify='center')
OwnAtmoEntry.pack(side=tk.LEFT)
def OwnAtmoEntryClear(event):
    OwnAtmoEntry.delete(0,"end")
OwnAtmoEntry.bind("<Button-1>", OwnAtmoEntryClear)
def OwnAtmoEntryFill(event):
    OwnAtmoEntry.insert(0,"0")
AtmoOwn = tk.IntVar()
def OwnAtmoEntryOut(x):
    AtmoOwn=OwnAtmoEntry.get()
    if is_num(AtmoOwn) == False:
        OwnAtmoEntryClear(x)
        OwnAtmoEntryFill(x)
OwnAtmoEntry.bind("<FocusOut>",OwnAtmoEntryOut)
LocationPredefFrame = tk.Frame(LocationFrame,
                               bd=3,
                               relief=tk.RIDGE)
LocationPredefFrame.pack(side=tk.LEFT)
AtmoOwnKmLabel = tk.Label(AtmoFrameOwn,
                          text="km").pack(side=tk.LEFT)
                        
#Display Window
DisplayFrame = tk.Frame(root,
                        bd=1,
                        padx=10,
                        pady=10,
                        relief=tk.SUNKEN)
DisplayFrame.pack(side=tk.BOTTOM)
DisplayText = tk.Text(DisplayFrame,
                      height=25,
                      width=130)
DisplayText.pack(side=tk.LEFT)
DisplayText.tag_configure('justcenter', justify='center')
scroll = tk.Scrollbar(DisplayFrame, command=DisplayText.yview)
scroll.pack(side=tk.RIGHT, fill=tk.Y)
DisplayText.configure(yscrollcommand=scroll.set)
quote = """
Welcome!
Purpose of this program is to calculate how relative atmosphere thickness changes during the day.
You will have to input a date, location and atmosphere thickness.
For location you can input your own or use one of the predefined locations.
As for the atmosphere thickness, you can choose multiple thicknesses (predefined and your own),
calculations will be done for each separetely. 
If you dont know which thickness to choose, click "Info" button for more details about Earth's atmosphere structure.
"""
DisplayText.insert(tk.END, quote, 'justcenter')

#Run Button
RunFrame = tk.Frame(root).pack()
buffer = io.StringIO()
def allya():
    buffer.seek(0)
    buffer.truncate(0)
    y=literal_eval(YearEntry.get())
    m=ComboMonths.current()+1
    d=literal_eval(ComboDays.get())
    RB = RBValue.get()
    AtmoOwn=OwnAtmoEntry.get()
    if RB == '1':
        phi=literal_eval(LatitudeEntry.get())
        lamb=literal_eval(LongitudeEntry.get())
        ltz=literal_eval(TimeZoneEntry.get())
    if RB == '2':
        CLV = ComboLocations.current()
        phi=literal_eval(Locations[CLV][2])
        lamb=literal_eval(Locations[CLV][4])
        ltz=literal_eval(Locations[CLV][6])
    atmo = []
    if var15.get() == 1:
        atmo.append(15)
    if var35.get() == 1:
        atmo.append(35)
    if var55.get() == 1:
        atmo.append(55)
    if var100.get() == 1:
        atmo.append(100)
    if CheckVar.get() == 1:
        atmo.append(literal_eval(AtmoOwn))

#Calculation Part
    def hms(x): #converts 0.0 to hh:mm:ss
        xh=x*24
        xm=(xh-int(xh))*60
        xs=(xm-int(xm))*60
        xSTRING=str("%2d"%int(xh))+':'+str("%02d"%int(xm))+':'+str("%02d"%int(xs))
        return xSTRING
        
    def jd2sza(zeroone,tstbase): #all from jd till solar zenital angle
        jd = julianday+tstbase
        jc = (jd-2451545)/36525 #Julian Century
        gmls = (280.46646+jc*(36000.76983+jc*0.0003032))%360 #Geometric Mean Longitute of the Sun (deg)
        gmas = 357.52911+jc*(35999.05029-0.0001537*jc) #Geometric Mean Anomaly of the Sun (deg)
        sqc = math.sin(math.radians(gmas))*(1.914602-jc*(0.004817+0.000014*jc))+math.sin(math.radians(2*gmas))*(0.019993-0.000101*jc)+math.sin(math.radians(3*gmas))*0.000289 #Sun's Equation of the Center
        stl = gmls+sqc #Sun's True Longitute (deg)
        sta = gmas+sqc #Sun's True Anomaly (deg)
        srv = (1.000001018*(1-ecc*ecc))/(1+ecc*math.cos(math.radians(sta))) #Sun's Radial Vector (AUs)
        scl = stl-0.00569-0.00478*math.sin(math.radians(125.04-1934.136*jc)) #Sun's Corrected Longitude (deg)
        moe = 23+(26+((21.447-jc*(46.815+jc*(0.00059-jc*0.001813))))/60)/60 #Mean Obliquity of the Ecliptic (deg)
        oc = moe+0.00256*math.cos(math.radians(125.04-1934.136*jc)) #Obliquity Correction (deg)
        sra = math.degrees(math.atan2(math.cos(math.radians(scl)),math.cos(math.radians(oc))*math.sin(math.radians(scl)))) #Sun's Right Ascension (deg)
        sd = math.degrees(math.asin(math.sin(math.radians(oc))*math.sin(math.radians(scl)))) #Suns' Declination (deg)
        varg = math.tan(math.radians(oc/2))*math.tan(math.radians(oc/2)) #var gamma
        eot = 4*math.degrees(varg*math.sin(2*math.radians(gmls))-2*ecc*math.sin(math.radians(gmas))+4*ecc*varg*math.sin(math.radians(gmas))*math.cos(2*math.radians(gmls))-0.5*varg*varg*math.sin(4*math.radians(gmls))-1.25*ecc*ecc*math.sin(2*math.radians(gmas)))#Equation of Time (minutes)
        if zeroone == 1: #one time just to get SN
            SN = (720-4*lamb-eot+ltz*60)/1440 #solar noon (lst)
            hasdivider = math.cos(math.radians(90.833))/(math.cos(math.radians(phi))*math.cos(math.radians(sd)))-math.tan(math.radians(phi))*math.tan(math.radians(sd))
            return SN, hasdivider
        if zeroone == 0: #for all other calculations
            tst = (tstbase*1440+eot+4*lamb-60*ltz)%1440 #True Solar Time (min)
            if tst/4 < 0: #Hour Angle (deg)
                ha = tst/4+180
            else:
                ha = tst/4-180
            sza = math.degrees(math.acos(math.sin(math.radians(phi))*math.sin(math.radians(sd))+math.cos(math.radians(phi))*math.cos(math.radians(sd))*math.cos(math.radians(ha)))) #Solar Zenith Angle (deg)
            return sza

    def alfabeta(): #calculate atmo thickness
        Sea = 90-sza #Solar Elevation Angle (deg) 
        if Sea < 0:
            Sea = 0
        alfa = math.degrees(math.asin(R/(R+atmo[a]) * math.cos(math.radians(Sea))))
        beta = 180 - 90 - alfa - Sea
        Thick = math.sin(math.radians(beta))/math.sin(math.radians(alfa)) * R
        return Thick, Sea

    #julian day
    A = math.floor(y/100)
    B = math.floor(A/4)
    C = 2-A+B
    E = math.floor(365.25*(y+4716))
    F = math.floor(30.6001*(m+1))
    julianday = C+d+E+F-1524.5-ltz/24
    ecc = 0.016704 #eccentricity of the earths orbit
    R = 6371
    SN, hasdivider = jd2sza(1,0)
    SNstr = hms(SN)
    if -1 < hasdivider < 1:
        has = math.degrees(math.acos(hasdivider)) #HA sunrise (deg)
        SRT = SN-has*4/1440 #Sunrise Time (lst)
        SST = SN+has*4/1440 #sunset time (lst)
        sld = 8*has #Sunlight Duration (minutes)
        SRTstr = hms(SRT)
        SSTstr = hms(SST)
        FromDawnTillNoon = int(24*(SN-SRT))
        FromNoonTillDusk = int(24*(SST-SN))
    else:
        SRT = SST = 0
        FromDawnTillNoon = FromNoonTillDusk = 11
        if phi > 0:
            if hasdivider*phi < 0:
                SRTstr = 'polar day'
                SSTstr = 'polar day'
            else:
                SRTstr = 'polar night'
                SNstr = 'polar night'
                SSTstr = 'polar night'
        else:
            if hasdivider*phi < 0:
                SRTstr = 'polar night'
                SNstr = 'polar night'
                SSTstr = 'polar night'
            else:
                SRTstr = 'polar day'
                SSTstr = 'polar day'
    ListOfThickPlusNoon = numpy.zeros((len(atmo),FromNoonTillDusk))
    SeaDusks = numpy.zeros((len(atmo),FromNoonTillDusk))
    ListOfThickMinusNoon = numpy.zeros((len(atmo),FromDawnTillNoon))
    SeaDawns = numpy.zeros((len(atmo),FromDawnTillNoon))
    ThickSN = [0]*len(atmo)
    ThickSRT = [0]*len(atmo)
    ThickSST = [0]*len(atmo)
    SeaSRT = [0]*len(atmo)
    SeaNoon = [0]*len(atmo)
    SeaSST = [0]*len(atmo)
    NoonPlusstr = numpy.empty((len(atmo),FromNoonTillDusk),dtype='object')
    DawnPlusstr = numpy.empty((len(atmo),FromDawnTillNoon),dtype='object')
    if type(ltz) == int:
        print('\nLocation:   \u03C6 =',"%+.2f" % (phi)+'\N{DEGREE SIGN}','  \u03BB =',"%+.2f" % (lamb)+'\N{DEGREE SIGN}','  Time Zone:',"%+d" % (ltz),file=buffer)
    if type(ltz) == float:
        print('\nLocation:   \u03C6 =',"%+.2f" % (phi)+'\N{DEGREE SIGN}','  \u03BB =',"%+.2f" % (lamb)+'\N{DEGREE SIGN}','  Time Zone:',"%+.1f" % (ltz),file=buffer)
    print("Day:       ",str(d)+'-'+str(m)+'-'+str(y),file=buffer)
    print("Sunrise:   ",SRTstr,file=buffer)
    print("Solar Noon:",SNstr,file=buffer)
    print("Sunset:    ",SSTstr,file=buffer)

    for a in range(len(atmo)): #each loop for each given atmosphere thickness
        #calc for solar noon (thinest thickness)
        sza = jd2sza(0,SN)    
        ThickSN[a], SeaNoon[a] = alfabeta()
        
        #calc for dawn (horizon - biggest thicnkess)
        if SRT > 0 or SRT < 0 or SRTstr == 'polar day':
            sza = jd2sza(0,SRT)
            ThickSRT[a], SeaSRT[a] = alfabeta()
            
        #calc for dusk (horizon - biggest thicnkess)
        if SST > 0 or SST < 0 or SRTstr == 'polar day':
            sza = jd2sza(0,SST)
            ThickSST[a], SeaSST[a] = alfabeta()
        
        for i in range(FromNoonTillDusk): #calculate everything from noon adding 1hour each loop till dusk
            NoonPlus = (SN*24+(i+1))/24
            NoonPlusstr[a][i] = hms(NoonPlus)
            sza = jd2sza(0,NoonPlus)
            ListOfThickPlusNoon[a][i], SeaDusks[a][i] = alfabeta()
                
        for i in range(FromDawnTillNoon): #calculate everything from dawn adding 1hour each loop till noon
            DawnPlus = (SN*24-FromDawnTillNoon+(i))/24
            DawnPlusstr[a][i] = hms(DawnPlus)
            sza = jd2sza(0,DawnPlus)
            ListOfThickMinusNoon[a][i], SeaDawns[a][i] = alfabeta()
                
        if SRT > 0 or SRT < 0 or SRTstr == 'polar day':
            print("\nAtmosphere thickness = "+str(atmo[a])+"km",file=buffer)
            print('  Time\t\tSolar Elevantion Angle[\N{DEGREE SIGN}]     Thickness[km]\tThickness(noon normalized)\tThickness(zenital normalized)',file=buffer)
            print(SRTstr,'\t\t\t',"%5.2f"%(SeaSRT[a]),'\t\t\t',"%6.2f"%(ThickSRT[a]),'\t\t   ',"%6.2f"%(ThickSRT[a]/ThickSN[a]),'\t\t\t\t',"%6.2f"%(ThickSRT[a]/atmo[a]),file=buffer)
            for i in range(FromDawnTillNoon):
                print(DawnPlusstr[a][i],'\t\t\t',"%5.2f"%(SeaDawns[a][i]),'\t\t\t',"%6.2f"%(ListOfThickMinusNoon[a][i]),'\t\t   ',"%6.2f"%(ListOfThickMinusNoon[a][i]/ThickSN[a]),'\t\t\t\t',"%6.2f"%(ListOfThickMinusNoon[a][i]/atmo[a]),file=buffer)
            print(SNstr,'\t\t\t',"%5.2f"%(SeaNoon[a]),'\t\t\t',"%6.2f"%(ThickSN[a]),'\t\t   ',"%6.2f"%(ThickSN[a]/ThickSN[a]),'\t\t\t\t',"%6.2f"%(ThickSN[a]/atmo[a]),file=buffer)
            for i in range(FromNoonTillDusk):
                print(NoonPlusstr[a][i],'\t\t\t',"%5.2f"%(SeaDusks[a][i]),'\t\t\t',"%6.2f"%(ListOfThickPlusNoon[a][i]),'\t\t   ',"%6.2f"%(ListOfThickPlusNoon[a][i]/ThickSN[a]),'\t\t\t\t',"%6.2f"%(ListOfThickPlusNoon[a][i]/atmo[a]),file=buffer)
            print(SSTstr,'\t\t\t',"%5.2f"%(SeaSST[a]),'\t\t\t',"%6.2f"%(ThickSST[a]),'\t\t   ',"%6.2f"%(ThickSST[a]/ThickSN[a]),'\t\t\t\t',"%6.2f"%(ThickSST[a]/atmo[a]),file=buffer)
    output = buffer.getvalue()
    DisplayText.insert(tk.END, output)
    DisplayText.yview_pickplace("end")
    
RunButton = tk.Button(RunFrame,
                      text="RUN",
                      width=10,
                      bd=3,
                      command=allya,
                      relief=tk.GROOVE).pack()
root.mainloop()
