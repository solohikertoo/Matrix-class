class Matrix:
    """Class to store and compute properties and subspaces of a two-dimensional matrix and 
    provide addition and multiplication operations.
    To create Matrix class: x = Matrix(a) where a is a list or tuple containing
    containing only lists or tuples of the same length. 
    Method are: 
    getRow(int), returns a Matrix that is the given row (0 based)
    getColumn(int), returns a Matrix that is the given column (0 based)
    makeTranspose(), returns a Matrix that is the transpose
    makeReduced(), returns a Matrix that is the row echelon form
    makeRange(), returns a Matrix in which columns span the range of the matrix
    makeNullSpace(), returns a Matrix in which columns span the nullspace of the matrix
    addMat(Matrix), returns Matrix, add the input to the matrix
    addScalar(float or int), returns Matrix, adds scalar to every element
    scalarMult(float or int), returns Matrix, multiplies every element by scalar
    leftMult(Matrix), returns a Matrix that is the result of multiply by input on the left
    rightMult(Matrix), returns a Matrix that is the result of multiplying by input on the right
    solveSys(Matrix), if input is a column vector, solves system or prints 'no solution', matrix
    does not have to be square, returns a Matrix that is a column vector containing solution,if an
    infinite number of solutions, the result is the base to which a combination of null space vectors
    can be added
    printMat(), prints the matrix"""
    
    EPS = 1e-6
    def __init__(self,data):
        if self.verify(data):
            self.data = [row[:] for row in data]
        else:
            if data:
                print("Invalid initialization")
            self.data = []
		
    def isEmpty(self):
        if not self.data:
            return True
        else:
            return False
            
    def getRow(self,rownum):
        """return given row as Matrix"""
        if not self.data:
            return Matrix([])
        return Matrix([self.data[rownum]])
        
    def getColumn(self,colnum):
        """return given column as Matrix"""
        if not self.data:
            return Matrix([])
        return Matrix([[self.data[i][colnum]] for i in range(len(self.data))])
        
    def makeTranspose(self):
        """return Matrix that is the transpose of the matrix"""
        if not self.data:
            return Matrix([])
        return Matrix([[row[i] for row in self.data] for i in range(len(self.data[0]))])

    def makeReduced(self):
        """return Matrix that is reduced row echelon form of matrix"""
        r, pivotList = self.getReducedWithPivots()
        return Matrix(r)

    def makeRange(self):
        """return Matrix whose columns span the range of the matrix"""
        if not self.data:
            return Matrix([])
        numRows = len(self.data)
        r,pivots = self.getReducedWithPivots()
        #extract columns corresponding to pivot positions
        rs = []
        pivotCols = [p[1] for p in pivots]
        for i in range(numRows):
            rangeRow = [self.data[i][j] for j in pivotCols]
            rs.append(rangeRow)
        return Matrix(rs)
        
    def makeNullSpace(self):
        """create Matrix whose columns span nullspace"""
        if not self.data:
            return Matrix([])
        numCols = len(self.data[0])
        reduced,pivots = self.getReducedWithPivots()
        #get column numbers corresponding to pivots and free variables
        pivCol = [value[1] for value in pivots]
        freeCol = [i for i in range(numCols) if not i in pivCol]
        nullspan = []
        #for each free variable, add list corresponding to null space vector
        for fi in freeCol:
            nullrow = [0 for i in range(numCols)]
            nullrow[fi] = 1
            for p in pivots:
                nullrow[p[1]] = -1*reduced[p[0]][fi]
            nullspan.append(nullrow)
        #return transpose so result is in form of column vectors 
        return Matrix(nullspan).makeTranspose()
        
    def addMat(self,B):
        """add input to matrix, return empty if invalid input"""
        if not self.data:
            return Matrix([])
        b=B.data
        if not self.verify(b):
            return Matrix([])
        numRows = len(self.data)
        numCols = len(self.data[0])
        if not (numCols == len(b[0]) and numRows == len(b)):
            return Matrix([])
        p = [[self.data[i][j]+b[i][j] for j in range(numCols)] for i in range(numRows)]
        return Matrix(p)
        
    def addScalar(self,n):
        """add a scalar to all elements"""
        if not self.data:
            return Matrix([])
        if not isinstance(n,(int,float)):
            return Matrix([])
        return Matrix([[value+n for value in row] for row in self.data])
        
    def scalarMult(self,n):
        """multiply all elements by scalar"""
        if not self.data:
            return Matrix([])
        if not isinstance(n,(int,float)):
            return Matrix([])
        return Matrix([[value*n for value in row] for row in self.data])       
        
    def rightMult(self,B):
        """multipy matrix on the right by input, return empty if invalid input"""
        if not self.data:
            return Matrix([])
        b=B.data
        if not self.verify(b):
            return Matrix([])
        numRows = len(self.data)
        numCols = len(self.data[0])
        if not numCols == len(b):
            return Matrix([])
        p = []
        for row in self.data:
            prodRow = [sum([row[i]*b[i][j] for i in range(numCols)]) for j in range(len(b[0]))]
            p.append(prodRow)
        return Matrix(p)
        
    def leftMult(self,B):
        """multipy matrix on the left by input, return empty if invalid input"""
        if not self.data:
            return Matrix([])
        b=B.data
        if not self.verify(b):
            return Matrix([])
        numRows = len(b)
        numCols = len(b[0])
        if not numCols == len(self.data):
            return Matrix([])
        p = []
        for row in b:
            prodRow = [sum([row[i]*self.data[i][j] for i in range(numCols)]) for j in range(len(self.data[0]))]
            p.append(prodRow)
        return Matrix(p)      
        
    def solveSys(self,B):
        """solve system of equations formed by matrix and input column vector for rhs
        return empty if invalid input or no solution. If infinite solutions,
        the returned vector can be added to a combination of nullspace vectors"""
        if not self.data:
            return Matrix([])
        numRows = len(self.data)
        numCols = len(self.data[0])
        #if not (numRows == numCols):
        #    return Matrix([])
        b=B.data
        if not self.verify(b):
            return Matrix([])
        if not len(b) == numRows:
            return Matrix([])
        if not all([len(row) == 1 for row in b]):
            return Matrix([]) 
        #augment matrix and get solution by row reduction
        augmented = []
        for i in range(numRows):
            row = self.data[i][:]
            row.append(b[i][0])
            augmented.append(row)
        augMatrix = Matrix(augmented)    
        r,pivots = augMatrix.getReducedWithPivots()
        #if pivot in rhs column then no solution
        pivotCols = [value[1] for value in pivots]
        if numCols in pivotCols:
            print("no solution")
            return Matrix([])
        #list of rhs after reduction
        rhs = [row[-1] for row in r]
        soln = [[0]] * numCols
        for p in pivots:
            soln[p[1]] = [rhs[p[0]]]          
        return Matrix(soln)
        
    def printMat(self):
        """print elements of matrix"""
        if self.isEmpty():
            print('Empty')
        else:
            for i in range(len(self.data)):
                print(self.data[i])
        print('\n')
        
        
    ################################
    #helpers
    
    def verify(self,data):
        """return True if the data used to initialize the matrix is in the proper form of a 
        list of lists, and that data contains ints or floats"""
        if not data:
            return False
        if not isinstance(data,(tuple,list)):
            return False
        ok = True
        for row in data:
            ok = ok and isinstance(row,(tuple,list))
        if ok:
            firstLen = len(data[0])
            for row in data:
                ok = ok and (len(row) == firstLen)
                for elem in row:
                    ok = ok and isinstance(elem,(int,float))
        return ok
   
    def getReducedWithPivots(self):
        """return reduced row echelon form and list of pivots"""
        if not self.data:
            return ([],[])
        #make copy of the matrix
        r = [row[:] for row in self.data]
        numRows = len(r)
        numCols = len(r[0])
        #forward substitution
        pivotColNum = 0
        pivotRowNum = 0
        pivotList = []
        while (pivotRowNum < numRows) and (pivotColNum < numCols):
            #row operations left to perform on current column
            currentColumn = True
            if abs(r[pivotRowNum][pivotColNum]) < Matrix.EPS:
                #if pivot is zero, swap with row with nonzero pivot and currentColumn is true,
                #or no more nonzero rows in this column so currentColumn is false
                currentColumn = self.swapNonZeroRow(r,pivotRowNum,pivotColNum)
            if currentColumn:
                #store pivot position (in reverse order for use in backward substitition
                pivotList.insert(0,(pivotRowNum,pivotColNum))
                pivot = r[pivotRowNum][pivotColNum]
                #normalize by pivot and perform row operations, when done
                r[pivotRowNum] = [r[pivotRowNum][i]/pivot for i in range(numCols)]
                for rowNum in range(pivotRowNum+1,numRows):
                    r[rowNum] = [r[rowNum][i]-r[rowNum][pivotColNum]*r[pivotRowNum][i] for i in range(numCols)]
                #move to next row only if there was a pivot for this row
                pivotRowNum = pivotRowNum+1
            #move to next column
            pivotColNum = pivotColNum + 1
        #backward substitution
        for pivotCoord in pivotList:
            pivotRowNum,pivotColNum = pivotCoord
            for rowNum in range(pivotRowNum-1,-1,-1):
                r[rowNum] = [r[rowNum][i]-r[rowNum][pivotColNum]*r[pivotRowNum][i] for i in range(numCols)]
        return (r,list(reversed(pivotList)))
        
    def swapNonZeroRow(self,r,pivotRowNum,pivotColNum):
        """if there is a row with a nonzero pivot swap with it and return true,
        otherwise no more nonzero rows in this column so return false"""
        numRows = len(r)
        rowNum = pivotRowNum+1
        while (rowNum < numRows) and (abs(r[rowNum][pivotColNum]) < Matrix.EPS):
            rowNum = rowNum +1
        if rowNum == numRows:
            return False
        else:
            tmpRow = r[pivotRowNum]
            r[pivotRowNum] = r[rowNum]
            r[rowNum] = tmpRow
            return True
                