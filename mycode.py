#had to use false as a str to get around null rules
f_f = "False"

#gets a resultant (its pretty much just Pythagorean Theorum)
def getResultant(vel1, vel2, np=False):
    #a^2 + b^2 = c^2
    answer = ((vel1**2)+(vel2**2))**0.5
    #These lines show up in every other equation np is noPrint and is used for the diagnostics
    if np == False:
        return print(answer)
    else:
        return answer

#this is the force equation just F=ma where: F=force, m=mass, a=acceleration
def getForce(f=f_f, m=f_f, a=f_f, np=False):
    #it checks to see if one of the args is set to "False" for substitution
    if f == f_f:
        answer = m*a
    elif m == f_f:
        answer = f/a
    else:
        answer = f/m
    if not np:
        return print(answer)
    else:
        return answer

#Kinematic equations are a set of equations with 5 potential vars being: vf=Velocity final, vi=Velocity initial,
#a=Acceleration, t=Time, d=Distance
#There are four diferent equations each one missing one of the variable and that is how they are name with KinematicNo<var>

#This is the kinematic equation without distance: vf = vi + a * t
def kinematicNoD(vf=f_f, vi=f_f, a=f_f, t=f_f, np=False):
    if vf == f_f:
        answer = vi + (a * t)
    elif vi == f_f:
        answer = vf - (a * t)
    elif a == f_f:
        answer = (vf - vi) / t
    elif t == f_f:
        #this part just makes sure that you don't use the wrong equation and try zero division
        if a == 0:
            answer = "That is zero division and therefore impossible"
        else:
            answer = (vf - vi) / a
    if not np:
        return print(answer)
    else:
        return answer

#kinematic equation without acceleration: d = (vi + vf) / 2 * t
def kinematicNoA(d=f_f, vi=f_f, vf=f_f, t=f_f, np=False):
    if d == f_f:
        answer = ((vi + vf) / 2) * t
    elif vi == f_f:
        answer = (2 * (d / t)) - vf
    elif vf == f_f:
        answer = (2 * (d / t)) - vi
    elif t == f_f:
        answer = d / ((vi + vf)/ 2)
    if not np:
        return print(answer)
    else:
        return answer

#kinematic equation without final velocity: d = vi * t + 1/2 * a * t^2
def kinematicNoVf(d=f_f, vi=f_f, t=f_f, a=f_f, np=False):
    if d == f_f:
        answer = (vi * t) + 0.5* a * (t **2)
    elif vi == f_f:
        answer = (d / t)-(0.5 * a * (t **2))
    elif t == f_f:
        #if the wrong numbers are fed in its very easy to get an imaginary number
        potential = ((-vi) + (((vi ** 2) - (4 * (0.5 * a) * (-d))) ** 0.5)) / (2 * (0.5 * a))
        #python imaginary numbers us "j" so if it has j in it it means it imaginary
        if "j" in str(potential):
            #this is an error message 
            answer = "You gave me the wrong number so now you must die (plz check your data)"
        elif potential < 0:
            #since the quadratic formula has two potential symobls (±) if it is less that zero the sign has to be changed
            answer = ((-vi) - (((vi **2 ) - (4 * (0.5 * a) * (-d))) ** 0.5)) / (2 * (0.5 * a))
        else:
            answer = potential
            answer
    elif a == f_f:
        answer = ((d / (vi * t)) / (0.5)) / (t **2)
    if not np:
        return print(answer)
    else:
        return answer

#kinematic equation without time: vf^2 = vi^2 + 2 * a * d
def kinematicNoT(vf=f_f, vi=f_f, a=f_f, d=f_f, np=False):
    if vf == f_f and vi != f_f and a != f_f and d != f_f:
        answer = ((vi **2 ) + 2 * a * d ) **0.5
    elif vi == f_f:
        answer = ((vf ** 2) - (2 * a * d)) ** 0.5
    elif a == f_f:
        answer = (((vf ** 2) - (vi ** 2)) / 2) / d
    elif d == f_f:
        answer = ((vf ** 2) - (vi ** 2)) / ( 2 * a)
    if not np:
        return print(answer)
    else:
        return answer

#this equation is specifically for filling in a projectile chart based on given information it just uses the kinematic equations but slightly adjusted
#variable are: Vx=Velocity on the x axis, dx=Distance on the x axis, t=Time, vfy=Velocity final on the y axis, dy=Distance on the y axis
def fillInChart(Vx=f_f, dx=f_f, t=f_f, vfy=f_f, dy=f_f):
    #velocity initial on the y axis
    viy = 0
    #acceleration on the y axis
    ay = -10
    if not dy == f_f and not dx == f_f:
        t = ((-viy) - (((viy ** 2) - (4 * (0.5 * ay) * (-dy))) ** 0.5)) / (2 * (0.5 * ay))
        vx = (dx / t)
        vfy = viy + (ay * t)
        answer = print("Your time is " +str(t)+ " your final and initial velocity for the x axis are " +str(vx)+ " and the final velocity for your y axis is " +str(vfy * -1))
    elif not Vx == f_f and not dx == f_f:
        t = (dx / ((Vx + Vx) / 2))
        dy = (viy * t) + 0.5 * ay * (t ** 2)
        vfy = ((viy ** 2) + 2 * ay * dy) ** 0.5
        answer = print("Your time is " +str(t)+ " Your distance on the y plane in " +str(dy)+ " Your velocity final on the y plane is " +str(vfy * -1))
    elif not t == f_f and not dx == f_f:
        Vx = (dx / t)
        dy = (viy * t) + 0.5 * ay * (t ** 2)
        vfy = ((viy ** 2) + 2 * ay * dy)** 0.5
        answer = print("Your velocity along the x axis is " +str(Vx)+ " You distance along the y axis is " +str(dy)+ " your velocity final along the y planes is " +str(vfy * -1) )
    elif not t == f_f and not Vx == f_f:
        dx = ((Vx + Vx) / 2) * t
        dy = (viy * t) + 0.5 * ay * (t ** 2)
        vfy = ((viy ** 2) + 2 * ay * dy) ** 0.5
        answer  = print("The distance traveled on the x plane is " +str(dx)+ " the distance traveled along the y plane is " +str(dy)+ " the final velocity along the y plane is " +str(vfy * -1))
    elif not Vx == f_f and not vfy == f_f:
        dy = ((vfy ** 2) - (viy ** 2)) / (2 * ay)
        t = dy / ((viy + vfy) / 2)
        dx = ((Vx + Vx) / 2) * t
        answer = print("The distance traveled along the y axis is " +str(dy)+ " the time taken is " +str((t**2)**0.5)+ " The distance traveled along the x axis is " +str(dx))
    elif not vfy == f_f and not dx == f_f:
        dy = ((vfy ** 2) - (viy ** 2)) / (2 * ay)
        t = dy / ((viy + vfy) / 2)
        Vx = (dx/t)
        answer = print("The distance traveled along the y axis is " +str(dy)+ " the time is " +str((t**2)**0.5)+ " the velocity along the x axis is " +str(Vx))
    elif not dy == f_f and not Vx == f_f:
        vfy = ((viy ** 2) + 2 * ay * dy) ** 0.5
        t = dy / ((viy + vfy) / 2)
        dx = ((Vx + Vx) / 2) *t
        answer = print("The final velocity along the y axis is " +str(vfy * -1)+ " the time is " +str((t**2)**0.5)+ " the distance along the x axis is " +str((dx **2)**0.5))
    return answer

#impulse equation: j = f * t
#j = impulse, f = force, t = time
def impulseEquation(j=f_f, f=f_f, t=f_f, np=False):
    if j == f_f:
        answer = f * t
    elif f == f_f:
        answer = j / t
    elif t == f_f:
        answer = j / f
    if not np:
        return print(answer)
    else:
        return answer

#momentum equation: p = m * v
#p=Momentum, m=Mass, vi=Velocity initial (typically zero but is adjustable), vf=Velocity final
def momentumEquation(p=f_f, m=f_f, vi=0, vf=f_f, np=False):
    if p == f_f:
        answer = m * (vf - vi)
    elif m == f_f:
        answer = p / (vf - vi)
    elif vf == f_f:
        answer =  p / m
    if not np:
        return print(answer)
    else:
        return answer

#this combines both of those: f * t = m * ∆v
#f = force, t=Time, m=Mass, ∆v = change in velocity
def impulseMomentumEquation(f=f_f, t=f_f, m=f_f, vi=0, vf=f_f, np=False):
    if f == f_f:
        answer = (m * (vf - vi)) / t
    elif t == f_f:
        answer = (m * (vf -vi)) / f
    elif m == f_f:
        answer = (f * t) / (vf - vi)
    elif vf == f_f:
        answer = ((f * t) / m) + vi
    if not np:
        return print(answer)
    else:
        return answer

#this one finds the arctan of an angle
#the equation is really complicated so you should probably look up the taylor series
#it is not fully acurate if you get within a 100th of 1 admitedly but it still works
def findTheta(opp, adj, np=False):
    #this is pi
    pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
    newvar = 0
    n = 0
    x = opp / adj
    nx = 1/x
    if x <= 1:
        for n in range(1000):
            var = ((-1)**n) * ((x**(2*n+1))/(2*n + 1))

            newvar += var
        newvar = newvar * (180/pi)
        if not np:
            print(newvar)
        else:
            return newvar
    else:
        for n in range(1000):
            var = ((-1)**n) * ((nx**(2*n+1))/(2*n + 1))
            newvar += var
        newvar = (newvar - (pi/2)) * -1
        newvar = newvar * (180/pi)
        if not np:
            print(newvar)
        else:
            return newvar

#conservation of momentum elastic: m1 * vi1 + m2 * vi2 = m1 * vf1 + m2 * vf2
#m1 = mass one, m2 = mass two, vi1=Velocity inital one, vi2=Velocity initial two, vf1=Velocity final one, vf2 = Velocity final two
def conservationOfMomentumE(m1, m2, v1i=f_f, v2i=f_f, v1f=f_f, v2f=f_f, np=False):
    if v2f == f_f:
        answer = (((m1 * v1i) + (m2 * v2i)) - (m1 * v1f)) / m2
    elif v1f == f_f:
        answer = (((m1 * v1i) + (m2 * v2i)) - (m2 * v2f)) / m1
    elif v1i == f_f:
        answer = (((m1 * v1f) + (m2 * v2f)) - (m2 * v2i)) / m1
    elif v2i == f_f:
        answer = (((m1 * v1f) + (m2 * v2f)) - (m1*  v1i)) / m2
    if not np:
        return print(answer)
    else:
        return answer

#conservation of momentum INelastic: m1 * vi1 + m2 * vi2 = vf(m1 + m2)
#its vf=Velocity final because in an inelastic equation they stick together and have the same final speed
def conservationOfMomentumI(m1, m2, v1i=f_f, v2i=f_f, vf=f_f, np=False):
    if vf == f_f:
        answer = ((m1 * v1i + (m2 * v2i))) / ((m1 + m2))
    elif v1i == f_f:
        answer = ((vf * (m1 + m2)) - (m2 * v2i)) / (m1)
    elif v2i == f_f:
        answer = ((vf * (m1 + m2)) - (m1 * v1i)) / (m2)
    if not np:
        return print(answer)
    else:
        return answer

#calculates work: W = F * d
#w=Work, f=Force, d=distance
def workEquation(w=f_f, f=f_f, d=f_f, np=False):
    if w == f_f:
        answer = f * d
    elif f == f_f:
        answer = w / d
    elif d == f_f:
        answer = w / f
    if not np:
        return print(answer)
    else:
        return answer

#finds kinetic energy: K = 1/2 * m * v^2
#ke = Kinetic Engergy, m=Mass, #velocity
def kineticEnergyEquation(ke=f_f, m=f_f, v=f_f, np=False):
    if ke == f_f:
        answer = (0.5 * m) * (v ** 2)
    elif m == f_f:
        answer = (ke / 0.5) / v ** 2
    elif v == f_f:
        answer = (ke / (0.5 * m)) ** 0.5
    if not np:
        return print(answer)
    else:
        return answer

#combines the previous two: f * d = (vi - vf) * (m * 0.5)
def workKineticEquation(f=f_f, d=f_f, m=f_f, vf=f_f, vi=0, np=False):
    if f == f_f:
        answer = (((vf - vi) ** 2) * (m * 0.5)) / d
    elif d == f_f:
        answer = (((vf - vi) ** 2) * (m * 0.5)) / f
    elif m == f_f:
        answer = ((f * d) / 0.5) / (vf - vi) ** 2
    elif vf == f_f:
        answer = (((f * d) / (0.5 * m)) - vi ** 2)  ** 0.5
    elif vi == f_f:
        answer = (((f * d) / (0.5 * m)) - vf**  2) ** 0.5
    if not np:
        return print(answer)
    else:
        return answer

#potential gavitational energy: ug = m * g * h
#ug=Potential gravitational energy, m=Mass, g=gravitational acceleration, h=height/distance
def potGravEquation(ug=f_f, m=f_f, g=f_f, h=f_f, np=False):
    if ug == f_f:
        answer = m * g * h
    elif m == f_f:
        answer = ug / (g * h)
    elif g == f_f:
        answer = ug / (m * h)
    elif h == f_f:
        answer = ug / (g * m)
    if not np:
        return print(answer)
    else:
        return answer

#potetial elastic energy: ue = 1/2 * k * x^2
#ue=Potential elastic energy, k=spring constant, x=length stretched/compressed
def potElasticEquation(ue=f_f, k=f_f, x=f_f, np=False):
    if ue == f_f:
        answer = (k * 0.5) * (x ** 2)
    elif k == f_f:
        answer = (ue / (x ** 2)) / 0.5
    elif x == f_f:
        answer = (ue / (k * 0.5)) ** 0.5
    if not np:
        return print(answer)
    else:
        return answer

#speed of circle equation: s = (2 * pi * d) / t
#s=Speed, d=distance, t=Time
def circularEquation(s=f_f, d=f_f, t=f_f, np=False):
    pi = 3.14
    if s == f_f:
        answer = (2 * pi * d) / t
    elif d == f_f:
        answer = (s * t) / (2 * pi)
    elif t == f_f:
        answer = (d * pi * 2) / s
    if not np:
        return print(answer)
    else:
        return answer

# TESTERS ===================================================================================================================================================
#all of these run diagnostics on the equations by testing for certian vars and ensuring that they are what they should be
#the smaller equations can be run to narrow down problems and the larger "runDiagnostics" function will tell you roughly where they are 
def testCircularEquation():
    s = circularEquation(d=0.25, t=1, np=True)
    d = circularEquation(s=1.57, t=1, np=True)
    t = circularEquation(s=0.785, d=0.25, np=True)
    if s == 1.57 and d == 0.25 and t == 2:
        diagnostics = True
    else:
        diagnostics = False
    return diagnostics

def testPotGravEquation():
    ug = potGravEquation(m=2, g=2, h=2, np=True)
    m = potGravEquation(ug=2, g=2, h=2, np=True)
    g = potGravEquation(ug=2, m=2, h=2, np=True)
    h = potGravEquation(ug=2, m=2, g=2, np=True)
    if ug == 8 and m == 0.5 and g == 0.5 and h == 0.5:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testPotElasticEquation():
    ue = potElasticEquation(k=2, x=2, np=True)
    k = potElasticEquation(ue=2, x=2, np=True)
    x = potElasticEquation(ue=4, k=2, np=True)
    if ue == 4 and k == 1 and x == 2:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testWorkKineticEquation():
    f = workKineticEquation(d=4.5, m=75, vf=6, np=True)
    d = workKineticEquation(f=4.5, m=75, vf=6, np=True)
    m = workKineticEquation(f=2, d=2, vf=2, np=True)
    vf = workKineticEquation(f=2, d=2, m=2, np=True)
    if f == 300 and d == 300 and m == 2 and vf == 2:
        diagnostic = True 
    else:
        diagnostic = False
    return diagnostic

def testKineticEnergyEquation():
    ke = kineticEnergyEquation(m=2, v=2, np=True)
    m = kineticEnergyEquation(ke=2, v=2, np=True)
    v = kineticEnergyEquation(ke=4, m=2, np=True)
    if ke == 4 and m == 1 and v == 2:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testWorkEquation():
    w = workEquation(f=2, d=2, np=True)
    f = workEquation(w=2, d=2, np=True)
    d = workEquation(w=2, f=2, np=True)
    if w == 4 and f == 1 and d == 1:
        diagnostic = True
    else:
        diagnostic = False 
    return diagnostic

def testConservationOfMomentumI():
    v1i = conservationOfMomentumI(m1=3, m2=6, v2i=5, vf=50, np=True)
    v2i = conservationOfMomentumI(m1=6, m2=3, vf=50, v1i=5, np=True)
    vf = conservationOfMomentumI(m1=1.5, m2=0.5, v1i=20, v2i=0, np=True)
    if v1i == 140.0 and v2i == 140.0 and vf == 15.0:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testConservationOfMomentumE():
    v2f = conservationOfMomentumE(m1=5, m2=7, v1i=4.5, v2i=-2, v1f=1.5, np=True)
    v1f = conservationOfMomentumE(m1=90, m2=10, v1i=0, v2i=0, v2f=-11, np=True)
    v1i = conservationOfMomentumE(m1=5, m2=6, v2i=3, v1f=6, v2f=9, np=True)
    v2i = conservationOfMomentumE(m1=5, m2=6, v1i=3, v1f=6, v2f=9, np=True)
    if v2f == 0.14285714285714285 and v1f == 1.2222222222222223 and v1i == 13.2 and v2i == 11.5:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testImpulseMomentumEquation():
    f = impulseMomentumEquation(t=2, m=2, vf=2, np=True)
    t = impulseMomentumEquation(f=2, m=2, vf=2, np=True)
    vf = impulseMomentumEquation(f=2, m=2, t=2, np=True)
    if f == 2 and t == 2 and vf == 2:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testFindTheta():
    test = findTheta(5, 1, True)
    if test == 78.69006752597979:
        diagnostic = True
    else:
        diagnostic = False

    return diagnostic

def testImpulseEquation():
    imp = impulseEquation(f=500, t=0.25, np=True)
    force = impulseEquation(t=3, j=900, np=True)
    time = impulseEquation(f=7500, j=18750, np=True)
    if imp == 125 and force == 300 and time == 2.5:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testMomentumEquation():
    moment = momentumEquation(m=1500, vf=35, np=True)
    mass = momentumEquation(p=17500, vf=5, np=True)
    vf = momentumEquation(p=97.5, m=65, np=True)
    if moment == 52500 and mass == 3500 and vf == 1.5:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testGetResultant():
    resultant = getResultant(3, 4, np=True)
    if resultant == 5:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testGetForce():
    force = getForce(m=2, a=2, np=True)
    mass = getForce(f=4, a=2, np=True)
    acc = getForce(f=4, m=2, np=True)
    if force == 4 and mass == 2 and acc == 2:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testKinematicNoD():
    vf = kinematicNoD(vi=2, a=2, t=2, np=True)
    vi = kinematicNoD(vf=2, a=2, t=2, np=True)
    a = kinematicNoD(vf=8, vi=2, t=2, np=True)
    t = kinematicNoD(vf=8, vi=2, a=2, np=True)
    if vf == 6 and vi == -2 and a == 3 and t == 3:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testKinematicNoA():
    d = kinematicNoA(vi=2, vf=2, t=2, np=True)
    vi = kinematicNoA(d=2, vf=2, t=2, np=True)
    vf = kinematicNoA(d=2, vi=2, t=2, np=True)
    t = kinematicNoA(d=2, vi=2, vf=2, np=True)
    if d == 4 and vi == 0 and vf == 0 and t == 1:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testKinematicNoVf():
    d = kinematicNoVf(vi=2, t=2, a=2, np=True)
    vi = kinematicNoVf(d=2, t=2, a=2, np=True)
    t = kinematicNoVf(d=2, vi=2, a=2, np=True)
    a = kinematicNoVf(d=2, vi =2, t=2, np=True)
    if d == 8 and vi == -3 and t == 0.7320508075688772 and a == 0.25:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

def testKinematicNoT():
    vf = kinematicNoT(vi=2, a=2, d=2, np=True)
    vi = kinematicNoT(vf=6, a=2, d=2, np=True)
    a = kinematicNoT(vf=6, vi=2, d=2, np=True)
    d = kinematicNoT(vf=6, vi=2, a=2, np=True)
    if vf == 3.4641016151377544 and vi == 5.291502622129181 and a == 8 and d == 8:
        diagnostic = True
    else:
        diagnostic = False
    return diagnostic

#this will tell you an overall result if anything isn't working
def runDiagnostics():
    diagnosticResults = ""
    if testPotGravEquation() and testPotElasticEquation() and testWorkKineticEquation() and testKineticEnergyEquation() and testWorkEquation() and testGetForce() and testGetResultant() and testKinematicNoA() and testKinematicNoD() and testKinematicNoT() and testKinematicNoVf() and testImpulseEquation() and testMomentumEquation() and testFindTheta() and testImpulseMomentumEquation() and testConservationOfMomentumE() and testConservationOfMomentumI() and testCircularEquation():
        diagnosticResults += "All good here"
    if not testGetForce():
        diagnosticResults += "getForce is malfuctioning "
    if not testGetResultant():
        diagnosticResults += "getResultant is malfuctioning "
    if not testKinematicNoA():
        diagnosticResults += "kinematicNoA is malfuctioning "
    if not testKinematicNoD():
        diagnosticResults += "kinematicNoD is malfuctioning "
    if not testKinematicNoT():
        diagnosticResults += "kinematicNoT is malfuctioning "
    if not testKinematicNoVf():
        diagnosticResults += "kinematicNoVf is malfuctioning "
    if not testImpulseEquation():
        diagnosticResults += "impulseEquation is malfuctioning "
    if not testMomentumEquation():
        diagnosticResults += "momentumEquation is malfuctioning "
    if not testFindTheta():
        diagnosticResults += "findTheta is malfuctioning "
    if not testImpulseMomentumEquation():
        diagnosticResults += "impulseMomentumEquation is malfuctioning "
    if not testConservationOfMomentumE():
        diagnosticResults += "conservationOfMomentumE is malfuctioning "
    if not testConservationOfMomentumI():
        diagnosticResults += "conservationOfMomentumI is malfuctioning "
    if not testWorkEquation():
        diagnosticResults += "workEquation is malfuctioning "
    if not testKineticEnergyEquation():
        diagnosticResults += "kineticEnergyEquation is malfuctioning "
    if not testWorkKineticEquation():
        diagnosticResults += "workKineticEquation is malfuctioning "
    if not testPotElasticEquation():
        diagnosticResults += "potElasticEquation is malfuctioning "
    if not testPotGravEquation():
        diagnosticResults += "potGravEquation is malfuctioning "
    if not testCircularEquation():
        diagnosticResults += "circularEquation is malfuctioning "

    return print(diagnosticResults)



# WHAT EQUATION =============================================================================================================================================

def num(num):
    try:
        return float(num)
    except ValueError:
        return num





#this allows it to be usable from the terminal and makes it more user friendly
def whatEquation():
    choice = input("Which equation do you need? ")
    if choice == "getForce":
        print("Please type False for the variable you are solving for")
        m = num(input("What is the mass? "))
        f = num(input("What is the force? "))
        a = num(input("what is the acceleration? "))
        getForce(m=m, f=f, a=a)
    elif choice == "kinematicNoT":
        print("PLEASE TYPE False FOR THE VARIABLE YOU ARE SOLVING FOR")
        vf = num(input("What is the final velocity? "))
        vi = num(input("What is the initial velocity? "))
        a = num(input("What is the acceleration? "))
        d = num(input("What is the distance? "))
        kinematicNoT(vf=vf, vi=vi, a=a, d=d)
    elif choice == "kinematicNoVf":
        print("PLEASE TYPE False FOR THE VARIABLE YOU ARE SOLVING FOR")
        d = num(input("What is the distance? "))
        vi = num(input("What is the initial velocity? "))
        t = num(input("What is the time? "))
        a = num(input("What is the acceleration? "))
        kinematicNoVf(d=d, vi=vi, t=t, a=a)
    elif choice == "kinematicNoA":
        print("PLEASE TYPE False FOR THE VARIABLE YOU ARE SOLVING FOR")
        d = num(input("What is the distance? "))
        vi = num(input("What is the initial velocity? "))
        vf = num(input("What is the final velocity? "))
        t = num(input("What is the time? "))
        kinematicNoA(d=d, vi=vi, vf=vf, t=t)
    elif choice == "kinematicNoD":
        print("PLEASE TYPE False FOR THE VARIABLE YOU ARE SOLVING FOR")
        vi = num(input("What is the initial velocity? "))
        vf = num(input("What is the final velocity? "))
        a = num(input("What is the acceleration? "))
        t = num(input("What is the time? "))
        kinematicNoD(vi=vi, vf=vf, a=a, t=t)
    elif choice == "getResultant":
        vel1 = num(input("What is the first number? "))
        vel2 = num(input("What is the second number? "))
        getResultant(vel1=vel1, vel2=vel2)
    elif choice == "fillInChart":
        vars = input("What variables do you have? ")
        if vars == "dy and dx":
            dx = num(input("What is your dx value? "))
            dy = num(input("What is your dy value? "))
            fillInChart(dy=dy, dx=dx)
        elif vars == "t and dx":
            t = num(input("What is you time? "))
            dx = num(input("What is your dx value? "))
            fillInChart(t=t, dx=dx)
        elif vars == "vfy and vx":
            vfy = num(input("What is your vfy value? " ))
            vx = num(input("What is your vx value? "))
            fillInChart(vfy=vfy, Vx=vx)
        elif vars == "vfy and dx":
            vfy = num(input("What is your vfy value? "))
            dx = num(input("What is your dx value? "))
            fillInChart(vfy=vfy, dx=dx)
        elif vars == "dy and vx":
            dy = num(input("What is your dy value? "))
            vx = num(input("What is your vx value? "))
            fillInChart(dy=dy, Vx=vx)
    elif choice == "run diagnostics":
        print("Running diagnostics")
        runDiagnostics()
    elif choice == "impulseEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        j = num(input("What is the impulse? "))
        f = num(input("What is the force? "))
        t = num(input("What is the time? "))
        impulseEquation(j=j, f=f, t=t)
    elif choice == "momentumEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        p = num(input("What is the momentum? "))
        m = num(input("What is the mass? "))
        vf = num(input("What is the final velocity? "))
        vi = num(input("What is the initial velocity (if you don't have one just put 0)? "))
        momentumEquation(p=p, m=m, vf=vf, vi=vi)
    elif choice == "impulseMomentumEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        f = num(input("What is you force? "))
        t = num(input("What is your time? "))
        m = num(input("What is your mass? "))
        vf = num(input("What is your final velocity? "))
        vi = num(input("What is your initial velocity (if you don't have this just put 0)? "))
        impulseMomentumEquation(f=f, t=t, m=m, vf=vf, vi=vi)
    elif choice == "findTheta":
        opp = num(input("What is the opposite side of the triangle? "))
        adj = num(input("What is the adjacent side of the triangle? "))
        findTheta(float(opp), float(adj))
    elif choice == "conservationOfMomentumE":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        m1 = num(input("What is the 1st mass? "))
        m2 = num(input("What is the 2nd mass? "))
        v1i = num(input("What is the 1st initial velocity? "))
        v2i = num(input("What is the 2nd initial velocity? "))
        v1f = num(input("What is the 1st final velocity? "))
        v2f = num(input("What is the 2nd final velocity? "))
        conservationOfMomentumE(m1=m1, m2=m2, v1i=v1i, v2i=v2i, v1f=v1f, v2f=v2f)
    elif choice == "conservationOfMomentumI":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        m1 = num(input("What is the 1st mass? "))
        m2 = num(input("What is the 2nd mass? "))
        v1i = num(input("What is the 1st initial velocity? "))
        v2i = num(input("What is the 2nd initial velocity? "))
        vf = num(input("What is the final velocity? "))
        conservationOfMomentumI(m1=m1, m2=m2, v1i=v1i, v2i=v2i, vf=vf)
    elif choice == "workEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        f = num(input("What is the force? "))
        d = num(input("What is the displacement? "))
        w = num(input("What is the work? "))
        workEquation(f=f, d=d, w=w)
    elif choice == "kineticEnergyEquation":
        print("PLEASE ENTER False FOR WHAY YOU ARE SOLVING FOR")
        ke = num(input("What is the Kinetic Energy? "))
        m = num(input("What is the mass? "))
        v = num(input("What is the veloctiy? "))
        kineticEnergyEquation(ke=ke, m=m, v=v)
    elif choice == "workKineticEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        f = num(input("What is the force? "))
        d = num(input("What is the displacement? "))
        m = num(input("What is the mass? "))
        vf = num(input("What is the final velocity? "))
        vi = num(input("What is the initial velocity (if you don't have this just put 0)? "))
        workKineticEquation(f=f, d=d, m=m, vf=vf, vi=vi)
    elif choice == "potGravEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        ug = num(input("What is your gravitational potential energy? "))
        m = num(input("What is your mass? "))
        g = num(input("What is the gravity? "))
        h = num(input("What is the height? "))
        potGravEquation(ug=ug, m=m, g=g, h=h)
    elif choice == "potElasticEquation":
        ue = num(input("What is your elastic potential energy? "))
        k = num(input("What is the k value? "))
        x = num(input("What is the x value? "))
        potElasticEquation(ue=ue, k=k, x=x)
    elif choice == "circularEquation":
        print("PLEASE ENTER False FOR WHAT YOU ARE SOLVING FOR")
        s = num(input("What is the speed? "))
        r = num(input("What is the radius? "))
        t = num(input("What is the period? "))
        circularEquation(s=s, d=r, t=t)
    else:
        print("Sorry that is not currently a valid equation. If you think it should be, feel free to make your own physics calculator from scratch")

#this is just called by default
whatEquation()

#hope you enjoyed the tour!








































