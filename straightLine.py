from Code.Python.robot_calc_functions import Magnitude, matmult, TransInv, Adjoint, FKinSpace, JacobianSpace, TransToRp, so3ToVec, VecToso3, MatrixLog3, RotInv
import numpy as np
from math import acos, tan, pi


# There was an error found in the source code for IKinSpace and IKinBody. The and condition in the if statement should be replaced with or
def IKinSpace(Slist, M, T, thetalist0, eomg, ev):
    maxiterations = 20
    success = False
    thf =[]#Final return variable
    thf.append(thetalist0)
    Vs = MatrixLog6(matmult(TransInv(FKinSpace(M, Slist, thetalist0)),T))
    wb  = Magnitude ([Vs[0],Vs[1],Vs[2]])
    vb = Magnitude ([Vs[3],Vs[4],Vs[5]])
    for i in range (maxiterations):
        if (wb>eomg or vb>ev):
            Jb = matmult(Adjoint(TransInv(FKinSpace(M, Slist, thetalist0))),JacobianSpace(Slist, thetalist0))
            thetalist0 = np.add(thetalist0, matmult(np.linalg.pinv(Jb),Vs))
            thf.append(thetalist0)
            Vs = MatrixLog6(matmult(TransInv(FKinSpace(M, Slist, thetalist0)),T))
            wb  = Magnitude ([Vs[0],Vs[1],Vs[2]])
            vb = Magnitude ([Vs[3],Vs[4],Vs[5]])
        else:
            success = True
            return (thf[len(thf)-1],success)
    return (thf[len(thf)-1],success)


def MatrixLog6(T):  # Takes a transformation matrix T SE(3)
    R, p = TransToRp(T)
    Rtrace = R[0][0] + R[1][1] + R[2][2]
    if (R == np.eye(3)).all():
        omg = [0, 0, 0]
        v = p
        theta = 1

    else:
        if (Rtrace == -1):
            theta = pi
            omg = MatrixLog3(R)
            G = (1 / theta) * np.eye(3) - 0.5 * np.asarray(VecToso3(omg)) + ((1 / theta) - (
            (1 / (tan(theta / 2.0))) / 2.0)) * (matmult(VecToso3(omg), VecToso3(omg)))
            v = np.dot(G, p)

        else:
            # Modified source code as I was running into issues with math domain errors from acos and a divide by zero error
            # with theta. I believe that these errors occur due to rounding and float calculation errors by the software. Because
            # of this, I am assuming that when theta is equal to zero, or the value in acos is > 1, it is in fact just a number close to
            # 0 or 1.
            if ((Rtrace - 1) / 2.0) > 1:
                Rtrace = 2.99999
            theta = acos((Rtrace - 1) / 2.0)
            if theta == 0:
                theta = .000001
            omg = so3ToVec((1 / (2 * np.sin(theta))) * (np.subtract(R, RotInv(R))))
            G = (1 / theta) * np.eye(3) - 0.5 * np.asarray(VecToso3(omg)) + ((1 / theta) - (
            (1 / (tan(theta / 2.0))) / 2.0)) * (matmult(VecToso3(omg), VecToso3(omg)))
            v = np.dot(G, p)

    return ([omg[0] * theta, omg[1] * theta, omg[2] * theta, v[0] * theta, v[1] * theta, v[2] * theta])


# createTraj is a function that will take in beginning and end coordinates to move the end effector of the UR5 6 joint robot from point to point.
# It will create a list of 101 points between these two points and use inverse kinematics to solve what joint angles are needed
# to get to the desired point in trajectory. This function uses the IKinSpace function; therefore the screws list used is with respect
# to the space frame.
#
def createTraj(begin_coords, end_coords):
    list_size = 101
    list_x = np.linspace(begin_coords[0], end_coords[0], list_size)
    list_y = np.linspace(begin_coords[1], end_coords[1], list_size)
    list_z = np.linspace(begin_coords[2], end_coords[2], list_size)

    anglelist = []
    L1 = .425
    L2 = .392
    W1 = .109
    W2 = .082
    H1 = .089
    H2 = .095
    M = [[-1, 0, 0, L1 + L2], [0, 0, 1, W1 + W2], [0, 1, 0, H1 - H2], [0, 0, 0, 1]]
    Slist = [[0, 0, 1, 0, 0, 0],[0, 1, 0, -H1, 0, 0],[0, 1, 0, -H1, 0, L1],[0, 1, 0, -H1, 0, L1+L2],[0, 0, -1, -W1, L1+L2, 0],[0, 1, 0, H2-H1, 0, L1+L2]]
    # Blist =[[0, 1, 0, W1+W2, 0, L1+L2], [0, 0, 1, H2, -L1-L2, 0], [0, 0, 1, H2, -L2, 0],[0, 0, 1, H2, 0, 0], [0, -1, 0, -W2, 0, 0],[0, 0, 1, 0, 0, 0]]
    eomg = 0.1
    ev = 0.01
    theta_guess = [.1, .1, .1, .1, .1, .1]
    for i in range(0,list_size):
        target_coords = (list_x[i],list_y[i],list_z[i])
        T_goal = [[1, 0, 0, target_coords[0]], [0, 1, 0, target_coords[1]], [0, 0, 1, target_coords[2]], [0, 0, 0, 1]]
        # result = IKinBody(Blist, M, T_goal, theta_guess, eomg, ev)
        result = IKinSpace(Slist, M, T_goal, theta_guess, eomg, ev)
        if result[1] == False:
            print "Error!"
            break
        else:
            print "Success!"
        anglelist.append(result[0])
        theta_guess = result[0]
    return anglelist

#writeToFile will go through the list of angle joints returned from createTraj function and write them into a csv file name angle.csv
def writeToFile(jointangles):
    f = open('angle.csv','w')
    for i in range(0, len(jointangles)):
        joint = jointangles[i]
        f.write(",".join(map(str, joint)) + "\n")
    f.close()


def main():
    begin_coords = (.7, 0, .1)
    end_coords = (0, -.3, .5)
    angleList = createTraj(begin_coords, end_coords)
    writeToFile(angleList)

if __name__=="__main__":
    main()