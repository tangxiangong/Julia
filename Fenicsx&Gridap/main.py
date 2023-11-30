from mpi4py import MPI
from dolfinx import mesh

nx, ny = 8, 8  # 剖分数量
domain = mesh.create_unit_square(MPI.COMM_WORLD, n, n, mesh.CellType.quadrilateral)