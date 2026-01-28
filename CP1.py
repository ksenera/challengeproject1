"""
PROJECT   : ChallengeProject1.py
PROGRAMMER: Kushika Senera
COURSE    : SFWRTECH 4MA3 - Numerical Linear Algebra and Numerical Optimization
INSTRUCTOR: Gagan Sidhu
"""

# A 𝑛 × 𝑛 Hilbert matrix has entries 𝐻ij = 1/ i + j -1, so form 
"""
|1   1/2 1/3 ...|
|1/2 1/3 1/4 ...|
|1/3 1/4 1/5 ...|
| :   :   :   ⋱|
"""
def generatorHb(n):
    # empty n x n matrix
    H = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(1.0/(i+j+1))
        H.append(row)

# using x, a n-vector with all entries equal to 1
    x = []
    for i in range(n):
        x.append(1.0)

# n-vector b=Hx 


# lower triangle system Lx = b
def forward_sub(L, b):
    # declare and initialize output vector x
    n = len(b)
    x = [0.0] * n
    # for j = 1 to n {loop over cols.} 
    # got length of b vector and put it in n to start loop
    for j in range(n):
        # xj = bj/Ljj {compute soln. component}
        x[j] = b[j] / L[j][j]

        # for i = j + 1 to n 
        # j starts at index 0 so j + 1 = 1 & n = 4 so from index 1 up to index 3 
        # range start to stop index n-1 => 4-1 = 3 
        for i in range(j+1,n):
            # bi = bi - Lijxj {update RHS}
            b[i] = b[i] - L[i][j] * x[j]
    return x

# upper triangle system Ux = b 
def back_sub(U, b):
    # declare and initialize output vector x
    n = len(b)
    x = [0.0] * n
    # for j = n to 1 {loop backwards over cols.} 
    for j in reversed(range(n)):
        # xj = bj/Ujj {compute soln. component}
        x[j] = b[j] / U[j][j]
        # for i = 1 to j - 1
        # range starts at 1 and stops at j-1
        # using range(stop) => j => j-1 
        for i in range(j):
            # bi = bi - uijxj {update RHS}
            b[i] = b[i] - U[i][j] * x[j]
    return x

# Gaussian Elimination using Partial Pivoting solving Ax = b
def gauss_elim(A_matrix, b_vector):
    # declare A system and b vector 
    n = len(A_matrix)
    # want to modify only local copies of A sys and b vec 
    A = [row[:] for row in A_matrix]
    b = b_vector[:]
    # A2 Note asks to print transformation at each step of elim + pivot
    print("Copy of A_matrix: ")
    for row in A:
        print(row)
    print(f"Copy of b_vector: {b}")

    # from ALGO 2.4 for pivoting steps
    # for k = 1 to n - 1 {loop over cols.}
    for k in range(n-1):
        # find index p |Apk| >= |Aik| for k <= i <= n {search for pivot in current col}
        p = k
        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[p][k]):
                p = i
        # if p != k then interchange rows k and p {interchange rows if necessary}
        if p != k:
            print(f"Swapping row {k} with row {p}")
            A[k], A[p] = A[p], A[k]
            b[k], b[p] = b[p], b[k]
            for row in A:
                print(row)
            print(f"Updated b: {b}\n")
        # in Algo 2.4 this line is line 2 of Algo 2.3 
        # if akk = 0 then stop {stop if pivot is zero}
        if A[k][k] == 0:  # if a_kk = 0 then stop
            print("Error: break if pivot is zero")
            break

        print(f"Elimination for Column {k}")
        # for j = k + 1 to n {compute multipliers for current col}
        for i in range(k+1, n):
            # mik = aik/akk 
            m_ik = A[i][k] / A[k][k]
            print(f"Multiplier m_{i}{k} = {m_ik:.2f}")

            A[i][k] = 0 

            # for j = k + 1 to n 
            for j in range(k+1, n):
                # for i = k + 1 to n -> nested for j in for i
                # {apply transformation to remaining submatrix}
                A[i][j] = A[i][j] - m_ik * A[k][j] 

            # solving Ly = b (Step 2 page 12 by updating b)
            b[i] = b[i] - m_ik * b[k]

        # printing as per assignment 
        for row in A:
            print([round(val, 2) for val in row])
        print(f"Updated b: {[round(val, 2) for val in b]}\n")

    x = back_sub(A, b)
    return x

if __name__ == "__main__":
    # ex 2.10 lower triangular matrix L
    L = [[4.0, 0.0, 0.0],
         [2.0, -2.0, 0.0],
         [1.0, 3.0, 4.0]]

    # RHS vector b
    b1 = [1.0, -2.0, 19.0]

    # Solve the system
    x_vector1 = forward_sub(L, b1)

    print(x_vector1)

    # ex 2.13 upper triangular matrix U
    U = [[1.0, 3.0, 4.0],
         [0.0, -2.0, 2.0],
         [0.0, 0.0, 4.0]]

    # RHS vector b
    b2 = [11.0, -2.0, 4.0]

    # Solve the system
    x_vector2 = back_sub(U, b2)

    print(x_vector2)

    #ans = gauss_elim(A_matrix, b_vector)

    #if ans:
        #print("FINAL SOLUTION x vector:")
        #print([round(val, 2) for val in ans])