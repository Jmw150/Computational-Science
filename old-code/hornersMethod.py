# not my code, just used an example
def horner(x0, *a): #does this even work?
    '''
        Horner's method is an algorithm to calculate a polynomial at
        f(x0) and f'(x0)

        x0 - The value to avaluate
        a - An array of the coefficients

        The degree is the polynomial is set equal to the number of coefficients
    '''
    n = len(a)

    y = a[0]
    z = a[0]
    for j in range(1, n - 1):
        y = x0 * y + a[j]
        z = x0 * z + y

    y = x0 * y + a[-1]

    print('P(x0) =', y)
    print('P\'(x0) =', z)

# not my code end


def hornersMethod(x_0, *a) : #broken?
    n = len(a)
    
    y = a[0]
    z = a[0]
    
    i = 0
    while (i < n) :
        y = x_0 * y + a[i]
        z = x_0 * z + y
        i += 1

    y = x_0 * y + a[-1]

    print("P(x_0)  = " + str(y))
    print("P\'(x_0) = " + str(z))

def poly_horner(A, x): #just an example, not my code
    p = A[-1]
    i = len(A) - 2
    while i >= 0:
        p = p * x + A[i]
        i -= 1
    return p


poly = ( 1, 0, 0) #x^2

x = 2

print (poly_horner(poly, x))
