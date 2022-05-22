import math


def func():
    z1 = float(input(" Enter z1: "))

    z2 = float(input(" Enter z2 (enter -3 to stop): "))

    if (z2 < 0):
        return 0

    s = (8 * pow(z1*z2, 1.5)) / pow(z1+z2, 3)                #Calculating overlap integrals
    
    h11 = 0.5*z1*z1 - 2*z1                                   #Calculating H matrix elements
    h22 = 0.5*z2*z2 - 2*z2
    h12 = pow(z1*z2, 1.5) * (4*z1*z2-8*z1-8*z2)/pow(z1+z2, 3)
    
                                                              #Calculating electron repulsion integral values
    r1111 = 5*z1/8
    r2222 = 5*z2/8
    d = pow(z1+z2, 4)
    r1122 = (pow(z1, 4)*z2 + 4*pow(z1, 3)*z2*z2 + z1*pow(z2, 4)+ 4*z1*z1*pow(z2, 3))/d
    r1212 = 20*pow(z1, 3)*pow(z2, 3)/pow(z1+z2,5)
    temp = (12*z1+8*z2)/pow(z1+z2, 2)+(9*z1+z2)/(2*z1*z1)
    r1112 = temp*16*pow(z1,4.5)*pow(z2,1.5)/pow(3*z1+z2,4)
    temp = (12*z2+8*z1)/pow(z1+z2,2)+(9*z2+z1)/(2*z2*z2)
    r1222 = temp*16*pow(z2,4.5)*pow(z1,1.5)/pow(3*z2+z1,4)

    k = float(input(" Guess for c1/c2 "))

    c2 = 1/math.sqrt(1+k*k+2*k*s)
    c1 = k*c2
    print(f" c1 = {c1}, c2 = {c2}")

    n = 0
    converged = 0

    while(converged == 0):
        p11 = 2*c1*c1                                       #Calculate probability density matrix elements
        p12 = 2*c1*c2
        p22 = 2*c2*c2
        
        f11 = h11+0.5*p11*r1111+p12*r1112+p22*(r1122-0.5*r1212)     #Calculate  matrix elements
        f12 = h12+0.5*p11*r1112+p12*(r1212*1.5-0.5*r1122)+0.5*p22*r1222
        f22 = h22+p11*(r1122-0.5*r1212)+p12*r1222+0.5*p22*r2222

        a = 1-s*s                               #Solve quadratic equation for the determinant
        b = 2*s*f12-f11-f22
        c = f11*f22-f12*f12
        rt = math.sqrt(b*b-4*a*c)

        e1 = (-b-rt)/(2*a)
        e2 = (-b+rt)/(2*a)

        e = e1

        if (e2 < e1):                    #Take lower value, solve roothian equation to get new coeffs d1 d2
            e=e2

        k = (e-f22)/(f12-s*e)
        d2 = 1/math.sqrt(1+k*k+2*k*s)
        d1 = k*d2

        n=n+1

        if(n > 500):
            print("Did not converge")
            return 0

        print(f" c1 = {c1}, c2 = {c2}, n = {n}")
        print(f" E = {e}")
        converged=1                         #converged flag
        if (abs(c2-d2) > 0.00001):          #not converged condition
            c1=d1
            c2=d2
            converged=0

        if (abs(c1-d1) > 0.00001):           #not converged condition
            c1=d1
            c2=d2
            converged=0



    if(converged == 1):
        print(f" Converged, N = {n}")
        ehf = e + 0.5*(p11*h11+2*p12*h12+p22*h22)   #Ehf obtained
        print(f" EHF = {ehf} Hartree")
        print(f" EHF = {27.2114*ehf} eV")


if __name__ == "__main__":
    func()
