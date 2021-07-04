# CPDINV

Software for estimating Curie-point depth using centroid method with both wavelet spectrum and Fourier spectrum. 
 
Curie-point depth is an important magnetic interface for studying the geological structure and thermal evolution of the crust and lithosphere, and it is common to estimate the Curie-point depth by spectral analysis of magnetic anomaly data. This FORTRAN software estimates Curie-point depth with the centroid method, using either Fourier transform or wavelet transform. It can be installed in Linux system with GFortran 4.8 and above.
This package consists of six functions, including projection, interpolation, power spectrum calculation, least squares fitting, Curie-point depth calculation.

Credits: Yihong Yin (15610046532@163.com) and Chun-Feng Li (cfli@zji.edu.cn) @ Zhejiang University

## **Installation**
### **Dependencies**
You will need GFORTRAN 4.8+. If your version is lower than it, you can refer to the update.sh file to update GFORTRAN version in system CentOS.
### **Installing using source code**
**1**. Download the repository from https://github.com/yinyih/CPDINV   
**2**. Go to your download directory,such as: `cd /CPDINV/code`  
**3**. Run the Makefile: `make install`  
**4**. Clear temporary files: `make clean`  
**5**. Add the program path under the code file directory to the PATH environment variable,such as:   
`echo “export  PATH=/home/CPDINV/code:$PATH”>>~/.bashrc`  
`source ~/.bashrc`  
**6**. Add documents,such as: `sudo cp /home/CPDINV/docs/*.1 /usr/share/man/man1`  
       If you want to know the detailed description of the command, use command such as `man Proj`

## **Uninstallation** 

If you want to remove or recompile the software, please use `make uninstall` command.  


## **References**

Cordell, L., Phillips, J.D., Godson, R.H., (1993). USGS potential-field geophysical software for PC and compatible microcomputers. Lead. Edge 12,290  

Kirby, J. F. (2005). Which wavelet best reproduces the Fourier power spectrum? Computers & Geosciences, 31(7), 846-864. doi:10.1016/j.cageo.2005.01.014  

Phillips, J.D., (1997). Potential-field geophysical software for the PC, version 2.2. U.S.Geological Survey, Open-File, Report, 97-725.  

Wang W Y.(2010). Minimum Curvature method and Fortran programming in data processing of potential field (in Chinese). Geological Publishing House.  
