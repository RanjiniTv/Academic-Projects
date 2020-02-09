
#Returns the total nuber of a particular type of cell in a grid
def SumCellType(iGrid, TypeOfCell):
    return sum(sum(iGrid == TypeOfCell))


# returns the position of the cell in the iRows*iCols matrix
def PointCoordinates(r, c, iRows, rankNum):
    i = rankNum * iRows + r - 1 
    return (i, c) 


#Distance from cancer cell at points to origin
def DistanceCal(points, origin):
    result = math.sqrt((points[0] - origin[0])**2 + (points[1] - origin[1])**2)
    return result


#Calculate Density of cancerous cells parameters: N prime and R
def DensityFormulaParam(grid, SubGridROWS):
    
    CancerCell = SumCellType(grid[1:SubGridROWS-1,:],'C')
    EffectorCell = SumCellType(grid[1:SubGridROWS-1,:],'E')
    DeadCell = SumCellType(grid[1:SubGridROWS-1,:],'D')
    
    NPrime = CancerCell + EffectorCell + DeadCell  
    RDistance = 0
    
    for r in range(1,SubGridROWS-1):
        for c in range(0,NumOfCols):
            if grid[r,c] == 'C':
                Coordinates = PointCoordinates(r, c, iRows, rank)	
                RDistance += DistanceCal(Coordinates, LatticeOrigin)
        return (c, NPrime, RDistance)	


# Calculate Density of cancerous cells: Rho(t) = N'(N prime)/R(squared)
def CancerDensityDev(nprime, r): 
    R = r / nprime
    return (nprime / (R**2))


# Determine the quadrant of various cell coordinates
def FindQuadrant(Coordinates, origin):
    
    if Coordinates[0] <= origin[0] and Coordinates[1] > origin[1]:
        return 'I'
    elif Coordinates[0] <= origin[0] and Coordinates[1] <= origin[1]:
        return 'II'
    elif Coordinates[0] > origin[0] and Coordinates[1] <= origin[1]:
        return 'III'
    elif Coordinates[0] > origin[0] and Coordinates[1] > origin[1]:
        return 'IV'


# calculates the growth of cells in each subgrid over a time 't'
def CellDevelopment(SubGrid, ProRate, DenseVal):    
    
    GridNew = numpy.copy(SubGrid)
    for r in range(0, SubROWS):
        for c in range(1, NumOfCols-1): 
            if SubGrid[r,c] == 'C':
                if random.random() < ProRate:
                    GridNew = mitosis(r, c, GridNew, LatticeOrigin, DenseVal) 
                elif random.random() < k2:
                    GridNew[r,c] = 'E'
            elif SubGrid[r,c] == 'E':
                if random.random() < k3:
                    GridNew[r,c] = 'D'
            elif SubGrid[r,c] == 'D':
                if random.random() < k4:
                    GridNew[r,c] = 'N'
    
    return GridNew





# Cell metastasis throughout the lattice
def mitosis(row, col, GridNew, LatticeOrigin, DenseVal):
    
    top = (row - 1, col)
    down = (row + 1, col)
    left = (row, col - 1) 
    right = (row, col + 1)
    AbnormalCells = ('C','E','D')
    
    InsideCells = {'I':[down,left],'II':[right,down],'III':[top,right],'IV':[left,top]}
    OutsideCells = {'I':[top,right],'II':[left,top],'III':[down,left],'IV':[right,down]}
    
    FirstInside = {'I':[down,left],'II':[right,down],'III':[right,right],'IV':[left,left]}
    FirstOutside = {'I':[right,right],'II':[left,left],'III':[down,left],'IV':[right,down]}
    
    LastInside = {'I':[left,left],'II':[right,right],'III':[top,right],'IV':[left,top]}
    LastOutside = {'I':[top,right],'II':[left,top],'III':[left,left],'IV':[right,right]}
     
    RandomVal = random.random()
    Coordinates = PointCoordinates(row, col, iRows, rank) 
    Quadrant = FindQuadrant(Coordinates, LatticeOrigin)
    
    if DenseVal:
        if row == 0:
            if RandomVal < 0.5 and GridNew[FirstOutside[Quadrant][0]] not in AbnormalCells:
                GridNew[FirstOutside[Quadrant][0]] = 'C'
            elif RandomVal >= 0.5 and GridNew[FirstOutside[Quadrant][1]] not in AbnormalCells:
                GridNew[FirstOutside[Quadrant][1]] = 'C' 
        elif row == (SubROWS-1):
            if RandomVal < 0.5 and GridNew[LastOutside[Quadrant][0]] not in AbnormalCells:
                GridNew[LastOutside[Quadrant][0]] = 'C'
            elif RandomVal >= 0.5 and GridNew[LastOutside[Quadrant][1]] not in AbnormalCells:
                GridNew[LastOutside[Quadrant][1]] = 'C'
        else:
            if RandomVal < 0.5 and GridNew[OutsideCells[Quadrant][0]] not in AbnormalCells:
                GridNew[OutsideCells[Quadrant][0]] = 'C'
            elif RandomVal >= 0.5 and GridNew[OutsideCells[Quadrant][1]] not in AbnormalCells:
                GridNew[OutsideCells[Quadrant][1]] = 'C'
    
    else:
        if row == 0:
            if RandomVal < 0.5 and GridNew[FirstInside[Quadrant][0]] not in AbnormalCells:
                GridNew[FirstInside[Quadrant][0]] = 'C'
            elif RandomVal >= 0.5 and GridNew[FirstInside[Quadrant][1]] not in AbnormalCells:
                GridNew[FirstInside[Quadrant][1]] = 'C'
        elif row == (SubROWS-1):
            if RandomVal < 0.5 and GridNew[LastInside[Quadrant][0]] not in AbnormalCells:
                GridNew[LastInside[Quadrant][0]] = 'C'
            elif RandomVal >= 0.5 and GridNew[LastInside[Quadrant][1]] not in AbnormalCells:
                GridNew[LastInside[Quadrant][1]] = 'C'
        else:
            if RandomVal < 0.5 and GridNew[InsideCells[Quadrant][0]] not in AbnormalCells:
                GridNew[InsideCells[Quadrant][0]] = 'C'
            elif RandomVal >= 0.5 and GridNew[InsideCells[Quadrant][1]] not in AbnormalCells:
                GridNew[InsideCells[Quadrant][1]] = 'C'
    return GridNew
        

#Transmission of data to higher ranks
def DataHighRank(SubGrid):
    
    comm.send(SubGrid[SubROWS-2,:],dest=rank+1)
    SubGrid[SubROWS-1,:]=comm.recv(source=rank+1)
    return 0

#Transmission of data to lower ranks
def DataLowRank(SubGrid):
    # send and receive rows with rank-1
    comm.send(SubGrid[1,:],dest=rank-1)
    SubGrid[0,:] = comm.recv(source=rank-1)
    return 0


if __name__ == '__main__':
    
    from mpi4py import MPI
    import time
    import numpy
    import math
    import random
    
    
    #MPI processing values
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    #Execution start time
    StartTime = MPI.Wtime()
    
    
    #Defining a square lattice of cells
    NumOfRows = 15
    NumOfCols = 15
    if size > NumOfRows:
        print("Number of processes is greater than row cells")
        exit()
        
    SubROWS = NumOfRows // size + 2
    iRows = size * (NumOfRows // size)
    iCols = NumOfCols
    LatticeOrigin = (iRows//2, iCols//2)
    
    
    # Inintialize a sub grid of normal cells and cancer cells
    SubGrid = numpy.empty([SubROWS, NumOfCols], dtype=object)
    SubGrid.fill('N')
    
    # At time t=0, cancer cells are present in the center of the grid
    if rank == (size//2):
        if size % 2 == 1:
            SubGrid[SubROWS//2, NumOfCols//2] = 'C'
            SubGrid[SubROWS//2-1, NumOfCols//2] = 'C'
            SubGrid[SubROWS//2+1, NumOfCols//2] = 'C'
            SubGrid[SubROWS//2, NumOfCols//2-1] = 'C'
            SubGrid[SubROWS//2, NumOfCols//2+1] = 'C'
        else:
            SubGrid[1, NumOfCols//2] = 'C'
            SubGrid[0, NumOfCols//2] = 'C'
            SubGrid[2, NumOfCols//2] = 'C'
            SubGrid[1, NumOfCols//2-1] = 'C'
            SubGrid[1, NumOfCols//2+1] = 'C'
    
    
    # Initial grid - Before cancer cell proliferation
    InitialGrid = comm.gather(SubGrid[1:SubROWS-1,:],root=0)
    if rank == 0:
        InitialGrid = numpy.vstack(InitialGrid)
        print("\n Initial Grid - Before cancer cell proliferation \n")
        print(InitialGrid[:])
      
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
        Param = DensityFormulaParam(SubGrid, SubROWS) 
        NumCancerCell = comm.gather(Param[0],root=0)
        nPrimeVal = comm.gather(Param[1],root=0)
        RVal = comm.gather(Param[2],root=0)
        DenseVal = False
        
        if rank == 0:
            NumCancerCell = sum(NumCancerCell)
            nPrimeVal = sum(nPrimeVal)
            RVal = sum(RVal)
            DenseVal =     (nPrimeVal, RVal) > RhoValue
        TotCancerCell = comm.bcast(NumCancerCell,root=0)
        DenseVal = comm.bcast(DenseVal,root=0)
        ProliferRate = k1 * (1 - TotCancerCell / PhiValue) 
        GridNew = CellDevelopment(SubGrid, ProliferRate, DenseVal) 
    	    
        if rank == 0:
            DataHighRank(GridNew)
        elif rank == size-1:
            DataLowRank(GridNew)
        else:
            DataHighRank(GridNew)
            DataLowRank(GridNew)
        SubGrid = numpy.copy(GridNew)
        
    
    # show final grid and number of each cell type present
    FinalGrid = comm.gather(SubGrid[1:SubROWS-1,:],root=0)
    if rank == 0:
        FinalGrid = numpy.vstack(FinalGrid)
        print("\n Final Grid - After cancer cell proliferation\n")
        print(FinalGrid[:])
        
        FinalNormalCells = SumCellType(FinalGrid,'N')
        FinalCancerCells = SumCellType(FinalGrid,'C')
        FinalEffectorCells = SumCellType(FinalGrid,'E')
        FinalDeadCells = SumCellType(FinalGrid,'D')
       
        print("\n Final Cells count")
        print("---------------------")
        print("Normal Cells: ",FinalNormalCells)
        print("Cancer cells: ",FinalCancerCells)
        print("Effector Cells: ",FinalEffectorCells)
        print("Dead Cells: ", FinalDeadCells)
        
        EndTime = MPI.Wtime()
        print("\n Total Parallel Execution Time:", EndTime - StartTime)