import matrix

x=matrix.Matrix([])
x.isEmpty()
x.printMat()

x=matrix.Matrix([[1,2,3],[5,6,7],[4,5,6]])
x.printMat()
x.getRow(0).printMat()
x.getColumn(1).printMat()
x.makeTranspose().printMat()
x.makeReduced().printMat()
x.makeRange().printMat()
x.makeNullSpace().printMat()

y=matrix.Matrix([[1,1,1],[2,2,2],[3,3,3]])
x.addMat(y).printMat()
x.rightMult(y).printMat()
x.leftMult(y).printMat()

x.addScalar(4).printMat()
x.scalarMult(2).printMat()

b=matrix.Matrix([[3],[11],[9]])
b.printMat()
x.solveSys(b).printMat()

x=matrix.Matrix([[1,0,2],[0,2,1],[0,3,3]])
x.printMat()
b=matrix.Matrix([[1],[1],[1]])
b.printMat()
x.solveSys(b).printMat()

x=matrix.Matrix([[1,0,1],[0,2,1],[0,3,3],[5,5,5]])
x.printMat()
b=matrix.Matrix([[1],[1],[1],[2]])
b.printMat()
x.solveSys(b).printMat()

x=matrix.Matrix([[1,0,2,4],[0,2,1,3],[0,1,1,1]])
x.printMat()
x.makeReduced().printMat()
x.makeNullSpace().printMat()
x.printMat()
b=matrix.Matrix([[1],[1],[1]])
b.printMat()
x.solveSys(b).printMat()
