import numpy as np
#Circle fitting itterative method, known as Taubin-Newton  
def Taubin_Newton (x_pos, y_pos):
    nPoints = len(x_pos)
    if (nPoints < 3):
        print("Too few points")
    else:
        x_centroid = np.mean(x_pos)
        y_centroid = np.mean(y_pos)
        Mxx = 0 
        Myy = 0
        Mxy = 0
        Mxz = 0
        Myz = 0
        Mzz = 0
        i = 0
        while (i < nPoints):
            Xi = x_pos[i] - x_centroid
            Yi = y_pos[i] - y_centroid
            Zi = Xi * Xi + Yi * Yi
            Mxy += Xi * Yi
            Mxx += Xi * Xi
            Myy += Yi * Yi
            Mxz += Xi * Zi
            Myz += Yi * Zi
            Mzz += Zi * Zi
            i = i+1

        Mxx /= nPoints
        Myy /= nPoints
        Mxy /= nPoints
        Mxz /= nPoints
        Myz /= nPoints
        Mzz /= nPoints

        Mz = Mxx + Myy
        Cov_xy = Mxx * Myy - Mxy * Mxy
        A3 = 4 * Mz
        A2 = -3 * Mz * Mz - Mzz
        A1 = Mzz * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz - Mz * Mz * Mz
        A0 = Mxz * Mxz * Myy + Myz * Myz * Mxx - Mzz * Cov_xy - 2 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy
        A22 = A2 + A2
        A33 = A3 + A3 + A3

        xnew = 0
        ynew = 1e+20
        epsilon = 1e-12
        iterMax = 20

        iter = 0
        while (iter < iterMax):
            yold = ynew
            ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3))
            if (np.abs(ynew) > np.abs(yold)):
                System.out.println("Newton-Taubin goes wrong direction: |ynew| > |yold|")
                xnew = 0
                break
            
            Dy = A1 + xnew * (A22 + xnew * A33)
            xold = xnew
            xnew = xold - ynew / Dy
            if (np.abs((xnew - xold) / xnew) < epsilon):
                break
            
            if (iter >= iterMax):
                System.out.println("Newton-Taubin will not converge")
                xnew = 0
            
            if (xnew < 0.):
                System.out.println("Newton-Taubin negative root: x = " + xnew)
                xnew = 0

            iter = iter + 1
                
        centreRadius = [0, 0, 0, 0]
        det = xnew * xnew - xnew * Mz + Cov_xy
        x = (Mxz * (Myy - xnew) - Myz * Mxy) / (det * 2)
        y = (Myz * (Mxx - xnew) - Mxz * Mxy) / (det * 2)
        centreRadius[0] = x + x_centroid
        centreRadius[1] = y + y_centroid
        centreRadius[3] = np.sqrt(centreRadius[0]**2 + centreRadius[1]**2)

        return centreRadius