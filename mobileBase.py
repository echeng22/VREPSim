from Code.Python.robot_calc_functions import FKinBody, TransInv, Adjoint, matmult, JacobianBody, MatrixLog6, Magnitude
import numpy as np
from math import pi, cos, sin
import matplotlib.pyplot as plt

def find_velocities(Blist, M0e, Tb0, F6, Vd, current_position):
    arm_pos = [current_position[3], current_position[4], current_position[5], current_position[6], current_position[7]]

    T0e = FKinBody(M0e, Blist, arm_pos)
    Jarm = JacobianBody(Blist, arm_pos)
    AdjVal = Adjoint(matmult(TransInv(T0e), TransInv(Tb0)))
    Jbase = matmult(AdjVal, F6)

    Je = [np.hstack((Jbase[0], Jarm[0])).tolist(),
          np.hstack((Jbase[1], Jarm[1])).tolist(),
          np.hstack((Jbase[2], Jarm[2])).tolist(),
          np.hstack((Jbase[3], Jarm[3])).tolist(),
          np.hstack((Jbase[4], Jarm[4])).tolist(),
          np.hstack((Jbase[5], Jarm[5])).tolist()]

    velocities = matmult(np.linalg.pinv(Je), Vd)
    return velocities


def find_current_EE_pos(Blist, M0e, Tb0, current_position):
    base_pos = [current_position[0], current_position[1], current_position[2]]
    arm_pos = [current_position[3], current_position[4], current_position[5], current_position[6], current_position[7]]

    Tsb = [[cos(base_pos[2]), -sin(base_pos[2]), 0, base_pos[0]],
           [sin(base_pos[2]), cos(base_pos[2]), 0, base_pos[1]],
           [0, 0, 1, .0963],
           [0, 0, 0, 1]]

    T0e = FKinBody(M0e, Blist, arm_pos)
    Tse = matmult(matmult(Tsb, Tb0), T0e)
    return Tse

def Ki_Integral(listXerr, currentErr):
    if len(listXerr) == 0:
        return np.asarray(currentErr) * .01
    else:
        sum_err = 0
        for i in range(len(listXerr)):
            sum_err = np.add(sum_err,.01 * np.asarray(listXerr[i]))
        return np.add(sum_err,currentErr)

def feedback(t, Xd, Xddot, X_current, listXerr):
    Kp = 10
    Ki = 5


    # print "test1"
    print "XD"
    print Xd
    print "XDDOT"
    print Xddot
    RRd = matmult(TransInv(Xd),Xddot)
    print RRd
    Vd = [RRd[2][1], RRd[0][2], RRd[1][0],RRd[0][3],RRd[1][3],RRd[2][3]]
    print Vd
    # print "test2"
    # print "X_CURRENT"
    # print X_current
    # print "XD"
    # print Xd
    # print "XDDOT"
    # print Xddot
    # print "VD"
    # print Vd
    print "----------------_ADJOINT----------------------"
    print Adjoint(matmult(TransInv(X_current), Xd))
    term1 = matmult(Adjoint(matmult(TransInv(X_current), Xd)),Vd)

    Xerror = MatrixLog6(matmult(TransInv(X_current),Xd))
    # print "-------------XERROR--------------"
    # print Xerror
    term2 = Kp*np.asarray(Xerror)
    term3 = Ki*Ki_Integral(listXerr, Xerror)
    # term3 = Ki*t*np.asarray(Xerror)
    # print "-------------TERM1--------------"
    # print term1
    # print "-------------TERM2--------------"
    # print term2
    # print "-------------TERM3--------------"
    # print term3
    sumTerms = np.add(np.add(term1, term2), term3)
    # print "-------------SUMTERMS--------------"
    # print sumTerms

    return (sumTerms, Xerror)


def X_traj(s):
    Xd = [[1, 0, 0, s],
          [0, sin(s*pi/2), cos(s*pi/2), 2*s + 1],
          [0, -cos(s*pi/2), sin(s*pi/2), .3 + .2*s],
          [0, 0, 0, 1]]
    return Xd

def Xdot_traj(t):
    sdot = sdot_time(t)
    s = s_time(t)
    sin_d = pi / 2 * cos(pi / 2 * s) * sdot
    cos_d = -pi / 2 * sin(pi / 2 * s) * sdot

    Xd = [[0, 0, 0, sdot],
          [0, sin_d, cos_d, 2*sdot],
          [0, -cos_d, sin_d, .2*sdot],
          [0, 0, 0, 0]]
    return Xd
    # return X_traj(sdot_time(t))

def sdot_time(t):
    return (6.0/25)*t - (6.0/125)*t**2

def s_time(t):
    return (3.0/25)*t**2 - (2.0/125)*t**3


def writeToFile(jointangles):
    f = open('traj.csv', 'w')
    for i in range(0, len(jointangles)):
        joint = jointangles[i]
        f.write(",".join(map(str, joint)) + "\n")
    f.close()

def main():
    r = .0475
    l = 0.47 / 2
    w = 0.3 / 2
    Blist = [[0, 0, 1, -.0330, 0, 0], [1, 0, 0, 0, -.5076, 0], [1, 0, 0, 0, -.3526, 0], [1, 0, 0, 0, -.2716, 0],
             [0, 0, 1, 0, 0, 0]]
    M0e = [[1, 0, 0, 0],
           [0, 1, 0, .0330],
           [0, 0, 1, .654],
           [0, 0, 0, 1]]
    Tb0 = [[1, 0, 0, 0],
           [0, 1, 0, .1662],
           [0, 0, 1, .0026],
           [0, 0, 0, 1]]
    F = (r / 4) * np.array([[-1 / (l + w), 1 / (l + w), 1 / (l + w), -1 / (l + w)],
                            [1.0, 1.0, 1.0, 1.0],
                            [-1.0, 1.0, -1.0, 1.0]])
    F6 = [[0, 0, 0, 0],
          [0, 0, 0, 0],
          [F[0][0], F[0][1], F[0][2], F[0][3]],
          [F[1][0], F[1][1], F[1][2], F[1][3]],
          [F[2][0], F[2][1], F[2][2], F[2][3]],
          [0, 0, 0, 0]]

    initpos = [0,0,0,0,0,-pi/2, pi/4, 0]
    list_pos = []
    list_pos.append(initpos)

    dt = .01
    t = 5
    listXerr = []
    for i in range(int(t/dt)):
        print i
        current_pos = list_pos[i]
        print current_pos
        time = i*dt
        s = s_time(time)
        print "-----------------------TIME VALUES------------------------------"
        print "-----------TIME---------------"
        print time
        print "-----------S---------------"
        print s
        X_current = find_current_EE_pos(Blist, M0e, Tb0, current_pos)
        Xd = X_traj(s)
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1Xd!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print Xd
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Current!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print X_current
        Xddot = Xdot_traj(time)

        V_t = feedback(time, Xd, Xddot, X_current, listXerr)
        listXerr.append(V_t[1])
        Vd = V_t[0]
        vel_list = find_velocities(Blist, M0e, Tb0, F6, Vd, current_pos)

        print "---------------VEL_LIST----------------"
        print vel_list

        wheel_vel = [vel_list[0], vel_list[1], vel_list[2], vel_list[3]]
        Vb = matmult(F6, wheel_vel)  # [0, 0, phi, x, y, z]
        print "---------------VB----------------"
        print Vb

        Vb3 = [Vb[2], Vb[3], Vb[4]]  # [phi, x, y]
        print "---------------VB3----------------"
        print Vb3
        if Vb3[0] != 0:
            wbz = Vb3[0]
            vbx = Vb3[1]
            vby = Vb3[2]
            Vb3[1] = (vbx * sin(wbz) + vby * (cos(wbz) - 1)) / wbz
            Vb3[2] = (vby * sin(wbz) + vbx * (1 - cos(wbz))) / wbz
        Tbs = [[1, 0, 0],
               [0, cos(current_pos[2]), -sin(current_pos[2])],
               [0, sin(current_pos[2]), cos(current_pos[2])]]
        Vb3 = matmult(Tbs, Vb3)
        update = [Vb3[1] * dt, Vb3[2] * dt, Vb3[0] * dt, vel_list[4] * dt, vel_list[5] * dt, vel_list[6] * dt,
                  vel_list[7] * dt, vel_list[8] * dt]

        next_pos = np.add(current_pos, update)
        list_pos.append(next_pos.tolist())
        print list_pos[len(list_pos) - 1]

    writeToFile(list_pos)

    # phi = []
    # x = []
    # y = []
    # Tau1 = []
    # Tau2 = []
    # Tau3 = []
    # Tau4 = []
    # Tau5 = []
    #
    # for i in range(0, len(list_pos)):
    #     x.append(list_pos[i][0])
    #     y.append(list_pos[i][1])
    #     phi.append(list_pos[i][2])
    #     Tau1.append(list_pos[i][3])
    #     Tau2.append(list_pos[i][4])
    #     Tau3.append(list_pos[i][5])
    #     Tau4.append(list_pos[i][6])
    #     Tau5.append(list_pos[i][7])

    phi_e = []
    x_e = []
    y_e = []
    Tau1_e = []
    Tau2_e = []
    Tau3_e = []

    for i in range(0, len(listXerr)):
        x_e.append(listXerr[i][0])
        y_e.append(listXerr[i][1])
        phi_e.append(listXerr[i][2])
        Tau1_e.append(listXerr[i][3])
        Tau2_e.append(listXerr[i][4])
        Tau3_e.append(listXerr[i][5])


    # timestamp = np.linspace(0, 5, 500)

    # plt.plot(timestamp, phi, label="phi")
    # plt.plot(timestamp, x, label="x")
    # plt.plot(timestamp, y, label="y")
    # plt.plot(timestamp, Tau1, label="Tau1")
    # plt.plot(timestamp, Tau2, label="Tau2")
    # plt.plot(timestamp, Tau3, label="Tau3")
    # plt.plot(timestamp, Tau4, label="Tau4")
    # plt.plot(timestamp, Tau5, label="Tau5")
    #
    # plt.ylim(-.2, .2)
    # plt.legend(loc='lower right')
    # plt.title("Plot of Configuration of Robot")
    # plt.show()

    plt.figure()
    timestamp = np.linspace(0, 5, 500)
    plt.plot(timestamp, phi_e, label="w1")
    plt.plot(timestamp, x_e, label="w2")
    plt.plot(timestamp, y_e, label="w3")
    plt.plot(timestamp, Tau1_e, label="v1")
    plt.plot(timestamp, Tau2_e, label="v2")
    plt.plot(timestamp, Tau3_e, label="v3")

    plt.ylim(-.2, .2)
    plt.legend(loc='lower right')
    plt.title("Plot of X_Error")
    plt.show()





if __name__=="__main__":
    main()
