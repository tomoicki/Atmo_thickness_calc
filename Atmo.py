import math, numpy
from ast import literal_eval

def is_num(x): #function which will check if input was a number
    try:
        float(x)
        return True
    except ValueError:
        return False
        
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
    
while True:
    y = input("Year: ")
    if is_num(y):
        break    
while True:
    m = input("Month (digit): ")
    if is_num(m) and 1 <= literal_eval(m) <= 12:
        break
while True:
    d = input("Day of the month: ")
    if is_num(d) and 1 <= literal_eval(d) <= 31:
        less31 = [4,6,9,11]
        if literal_eval(m) in less31:
            if literal_eval(d) > 30:
                print("The month you choose does not have that many days.")
                continue
        if literal_eval(m) == 2:
            if literal_eval(d) > 28:
                print("February does not have that many days.")
                continue
        break
while True:
    yesno = input('Latitude as for Wrocław (+51.1\N{DEGREE SIGN})? (enter y or n): ')
    if yesno == 'y':
        phi = '51.1'
        break
    if yesno == 'n':
        while True:
            phi = input('Latitude (+N, -S) (Decimal Degrees Format i.e. 0.00\N{DEGREE SIGN} not 0\N{DEGREE SIGN}'"0'"'0"): ')
            if is_num(phi) and -90 <= literal_eval(phi) <= 90:
                break
        break
while True:
    yesno = input('Longitute as for Wrocław (+17.03\N{DEGREE SIGN})? (enter y or n): ')
    if yesno == 'y':
        lamb = '17.03'
        break
    if yesno == 'n':
        while True:
            lamb = input('Longitude (+E, -W) (Decimal Degrees Format i.e. 0.00\N{DEGREE SIGN} not 0\N{DEGREE SIGN}'"0'"'0"): ')
            if is_num(lamb) and -180 <= literal_eval(lamb) <= 180:
                break
        break
ltz = round(literal_eval(lamb)/15)
while True:
    yesno = input("Time Zone "+str("%+d"%ltz)+"? (enter y or n): ")
    if yesno == 'y':
        ltz = str(ltz)
        break
    if yesno == 'n':
        while True:
            ltz = input("Time Zone (+E, -W)): ")
            if is_num(ltz) and -12 <= literal_eval(ltz) <= 12:
                break
        break
print("Calculations will be done assuming thickness of the Earth's atmosphere = 10km, 35km and 50km.")
while True:
    yesno = input("Want to add your own atmosphere thickness? (enter y or n): ")
    if yesno == 'n':
        atmo = [15,35,55]
        break
    if yesno == 'y':
        while True:
            r = input("Atmosphere thickness (km): ")
            if is_num(r):
                atmo = [10,35,50,literal_eval(r)]
                break
        break
phi = literal_eval(phi)
lamb = literal_eval(lamb)
ltz = literal_eval(ltz)
d = literal_eval(d)
m = literal_eval(m)
y = literal_eval(y)
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
print('\nLocation:   \u03C6 =',"%+.2f" % (phi)+'\N{DEGREE SIGN}','  \u03BB =',"%+.2f" % (lamb)+'\N{DEGREE SIGN}','  Time Zone:',"%+d" % (ltz))
print("Day:       ",str(d)+'-'+str(m)+'-'+str(y))
print("Sunrise:   ",SRTstr)
print("Solar Noon:",SNstr)
print("Sunset:    ",SSTstr)

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
        print("\nAtmosphere thickness = "+str(atmo[a])+"km")
        print('  Time\t\tSolar Elevantion Angle[\N{DEGREE SIGN}]     Thickness[km]\tThickness(noon normalized)\tThickness(zenital normalized)')
        print(SRTstr,'\t\t',"%5.2f"%(SeaSRT[a]),'\t\t\t',"%6.2f"%(ThickSRT[a]),'\t\t',"%6.2f"%(ThickSRT[a]/ThickSN[a]),'\t\t\t',"%6.2f"%(ThickSRT[a]/atmo[a]))
        for i in range(FromDawnTillNoon):
            print(DawnPlusstr[a][i],'\t\t',"%5.2f"%(SeaDawns[a][i]),'\t\t\t',"%6.2f"%(ListOfThickMinusNoon[a][i]),'\t\t',"%6.2f"%(ListOfThickMinusNoon[a][i]/ThickSN[a]),'\t\t\t',"%6.2f"%(ListOfThickMinusNoon[a][i]/atmo[a]))
        print(SNstr,'\t\t',"%5.2f"%(SeaNoon[a]),'\t\t\t',"%6.2f"%(ThickSN[a]),'\t\t',"%6.2f"%(ThickSN[a]/ThickSN[a]),'\t\t\t',"%6.2f"%(ThickSN[a]/atmo[a]))
        for i in range(FromNoonTillDusk):
            print(NoonPlusstr[a][i],'\t\t',"%5.2f"%(SeaDusks[a][i]),'\t\t\t',"%6.2f"%(ListOfThickPlusNoon[a][i]),'\t\t',"%6.2f"%(ListOfThickPlusNoon[a][i]/ThickSN[a]),'\t\t\t',"%6.2f"%(ListOfThickPlusNoon[a][i]/atmo[a]))
        print(SSTstr,'\t\t',"%5.2f"%(SeaSST[a]),'\t\t\t',"%6.2f"%(ThickSST[a]),'\t\t',"%6.2f"%(ThickSST[a]/ThickSN[a]),'\t\t\t',"%6.2f"%(ThickSST[a]/atmo[a]))