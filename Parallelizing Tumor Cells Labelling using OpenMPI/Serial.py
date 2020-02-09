
def FindQuadrant(row, col, LatticeOrigin):
    # get quadrant of coordinates relative to origin
    if row <= LatticeOrigin[0] and col > LatticeOrigin[1]:
        return 'I'
    elif row <= LatticeOrigin[0] and col <= LatticeOrigin[1]:
        return 'II'
    elif row > LatticeOrigin[0] and col <= LatticeOrigin[1]:
        return 'III'
    elif row > LatticeOrigin[0] and col > LatticeOrigin[1]:
        return 'IV'

#Distance from cancer cell at (row,col) to origin
def DistanceCal(row, col):
    result = math.sqrt((row - LatticeOrigin[0])**2 + (col - LatticeOrigin[1])**2)
    return result


#Returns the total nuber of a particular type of cell in a grid
def TotalCellType(TypeOfCell):
    return sum(sum(InitialGrid==TypeOfCell))


#Density of cancerous cells: Rho(t) = N'(N prime)/R(squared)
def CancerDensityGrowth(NumOfRows, NumOfCols):
    CancerCell = TotalCellType('C')
    DeadCell = TotalCellType('D')
    EffectorCell = TotalCellType('E')
    NPrime = CancerCell + DeadCell + EffectorCell
    RDistance = 0
    
    for r in range(NumOfRows):
        for c in range(NumOfCols):
            if InitialGrid[r, c] == 'C':
                RDistance += DistanceCal(r, c)
        RDistance = RDistance / NPrime
    return NPrime / RDistance**2
            

# Proliferation of cancer cells
def ProliferationRate(k, Nc, Phi):
    K1Prime = k * (1 - Nc / Phi)
    return K1Prime


# Cell metastasis throughout the lattice
def CancerMetastasis(RandomVal, DensityGrowth, Quadrant):
    
    top = (row - 1, col)
    down = (row + 1, col)
    left = (row, col - 1) 
    right = (row, col + 1)
    AbnormalCells = ('E','D')
    
    InsideCells = {'I':[down,left],'II':[right,down],'III':[top,right],'IV':[left,top]}
    OutsideCells = {'I':[top,right],'II':[left,top],'III':[down,left],'IV':[right,down]}
     
    
    if DensityGrowth:
        if RandomVal < 0.5 and ChangeGrid[OutsideCells[Quadrant][0]] not in AbnormalCells:
            ChangeGrid[OutsideCells[Quadrant][0]] = 'C'
        elif RandomVal >= 0.5 and ChangeGrid[OutsideCells[Quadrant][1]] not in AbnormalCells:
            ChangeGrid[OutsideCells[Quadrant][1]] = 'C'
    else:
        if RandomVal < 0.5 and ChangeGrid[InsideCells[Quadrant][0]] not in AbnormalCells:
            ChangeGrid[InsideCells[Quadrant][0]] = 'C' 
        elif RandomVal >= 0.5 and ChangeGrid[InsideCells[Quadrant][1]] not in AbnormalCells: 
            ChangeGrid[InsideCells[Quadrant][1]] = 'C'
            

if __name__ == '__main__':

    import numpy
    import random
    import math
    import time

    #Execution start time
    StartTime = time.time()
    
    
    #Defining a square lattice of cells
    NumOfRows = 15
    NumOfCols = 15
    LatticeOrigin = (NumOfRows//2, NumOfCols//2)
    
    #Inintialize a grid of normal cells and cancer cells
    InitialGrid = numpy.empty([NumOfRows,NumOfCols],dtype=object)
    InitialGrid.fill('N')
    
    #At time, t=0 , there are 5 cancer cells in center of grid
    InitialGrid[NumOfRows//2, NumOfCols//2] = 'C'
    InitialGrid[NumOfRows//2 + 1, NumOfCols//2] = 'C'
    InitialGrid[NumOfRows//2 - 1, NumOfCols//2] = 'C'
    InitialGrid[NumOfRows//2, NumOfCols//2 - 1] = 'C'
    InitialGrid[NumOfRows//2, NumOfCols//2 + 1] = 'C'
    
    print("\n Initial Grid - Before cancer cell proliferation \n")
    print(InitialGrid)
    
    ChangeGrid = numpy.copy(InitialGrid)
    
    
    #Constants in the proliferation or dissolution of cancer cells
    k1 = 0.7
    k2 = 0.2
    k3 = 0.3
    k4 = 0.3
    PhiValue = NumOfRows * NumOfCols    #Total number of cancer cells reaches a maximum
    
    #Density of cancer cells over a time 't'
    TimeLapse = 10
    RhoValue = 3.85                     #Density of cancer cells
    
   
    # Cancer cell growth over a time period 't'
    for t in range(TimeLapse):
        DensityGrowth = (CancerDensityGrowth(NumOfRows, NumOfCols) > RhoValue) 
        for row in range(1, NumOfRows):
            for col in range(1, NumOfCols):
                Quadrant = FindQuadrant(row, col, LatticeOrigin)
                if InitialGrid[row, col] == 'C':
                    if random.random() < ProliferationRate(k1, TotalCellType('C'), PhiValue):
                        RandomVal = random.random()
                        CancerMetastasis(RandomVal, DensityGrowth, Quadrant) 
                    elif random.random() < k2:
                        ChangeGrid[row, col] = 'E'
                elif InitialGrid[row, col] == 'E':
                    if random.random() < k3:
                        ChangeGrid[row, col] = 'D'
                elif InitialGrid[row, col] == 'D':
                    if random.random() < k4:
                        ChangeGrid[row, col] = 'N'
        InitialGrid = numpy.copy(ChangeGrid)
        
        
# Final grid after cell proliferation
print("\n Final Grid - After cancer cell proliferation\n")
print(InitialGrid)
FinalNormalCells = TotalCellType('N')
FinalCancerCells = TotalCellType('C')
FinalEffectorCells = TotalCellType('E')
FinalDeadCells = TotalCellType('D')

print("\n Final Cells count")
print("---------------------")
print("Normal Cells: ",FinalNormalCells)
print("Cancer cells: ",FinalCancerCells)
print("Effector Cells: ",FinalEffectorCells)
print("Dead Cells: ", FinalDeadCells)

EndTime = time.time()
print("\nTotal Serial Execution Time:", EndTime - StartTime)