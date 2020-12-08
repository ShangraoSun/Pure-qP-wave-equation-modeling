# Pure-qP-wave-equation-modeling
Finite Difference (FD) method  is used to solve the pure qP-wave equation; the methods of efficiently solving the implicit FD stencil (Poisson equation) involved here can be discussed

Note that, before run the code - 'weakaniso_qP_forward.m', one should creat the file 'qP_Datas' for saving files created by this code.

# Codes Description
The main program named 'weakaniso_qP_forward.m', is about the pure qP-wave euqation for TTI media. 
- Large Linear Sparse Matrix Equation
I show how to creat Sparse Matrix in the main program (for solving implicit FD eqaution), and place the 2nd-orders of spatial FD. Also, I creat bianry and matlab files which are saved in 'qP_Datas' file, to record the models used here. 
- Absorbing Boundary Conditions (ABCs)
I choose the exponential ABCs proposed by Cerjan (1985) to absorb the boundary reflection. By using this ABCs, calculation efficiency will drop further.
- 'extendmodel' program
This program is aim to extend models used in our main program.

# Discussion
Although solving implicit FD equation directly is not a good choice, I still want to find a suitable solving method for it. After having a research for solution of Poisson equation , I will try to use 'Fast Poisson Solver' (propoed by Stang, 2007) to solve the Poisson equation. Any good method is welcome.

Any quastion, touch me: 13152175221@163.com
