from Code.Python.robot_calc_functions import ForwardDynamics, EulerStep
#from robot_calc_functions import ForwardDynamics, EulerStep
import numpy as np
import matplotlib.pyplot as plt

# Needed to edit this function as the original function would add all steps solved by Euler to thetamat, instead of the
# final resolution step of Euler's Integration. Shifted the thetamat.append and dthetamat.append back one tab to be inside
# the first for loop, not in the 2nd for loop.
def ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes):
    if intRes < 1 or isinstance(intRes, float):
        print "Integration resolution must be an integer value greater than zero."
        return
    if Ftipmat == 0:
        NewFtipmat = [[0,0,0,0,0,0]]*len(taumat)
    else:
        NewFtipmat = Ftipmat
    thetamat = []
    thetamat.append(thetalist)
    dthetamat = []
    dthetamat.append(dthetalist)
    for i in range(len(taumat)):
        for j in range(intRes):
            ddthetalist = ForwardDynamics(thetalist, dthetalist, taumat[i], g, NewFtipmat[i], Mlist, Glist, Slist)
            thetalist,dthetalist = EulerStep(thetalist,dthetalist,ddthetalist,(dt/intRes))
        thetamat.append(thetalist)
        dthetamat.append(dthetalist)
    return thetamat, dthetamat

def fallingSim():
    thetalist = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    dthetalist = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Steiner's based on what was given on URDF XYZ change
    G1 = [[0.010267495893, 0, 0, 0 ,0, 0], [0, 0.010267495893, 0, 0 ,0, 0], [0, 0, .00666, 0 ,0, 0], [0, 0, 0, 3.7 ,0, 0], [0, 0, 0, 0 ,3.7, 0], [0, 0, 0, 0, 0, 3.7]]
    G2 = [[0.22689067591 + 8.393*.28**2, 0, 0, 0, 0, 0], [0, 0.22689067591 + 8.393*.28**2, 0, 0, 0, 0], [0, 0, .0151074, 0, 0, 0], [0, 0, 0, 8.393, 0, 0], [0, 0, 0, 0, 8.393, 0], [0, 0, 0, 0, 0, 8.393]]
    G3 = [[0.049443313556 + 2.275*.25**2, 0, 0, 0, 0, 0], [0, 0.049443313556  + 2.275*.25**2, 0, 0, 0, 0], [0, 0, 0.004095, 0, 0, 0], [0, 0, 0, 2.275, 0, 0], [0, 0, 0, 0, 2.275, 0], [0, 0, 0, 0, 0, 2.275]]
    G4 = [[0.111172755531, 0, 0, 0, 0, 0], [0, 0.111172755531, 0, 0, 0, 0], [0, 0, 0.21942, 0, 0, 0], [0, 0, 0, 1.219, 0, 0], [0, 0, 0, 0, 1.219, 0], [0, 0, 0, 0, 0, 1.219]]
    G5 = [[0.111172755531, 0, 0, 0, 0, 0], [0, 0.111172755531, 0, 0, 0, 0], [0, 0, 0.21942, 0, 0, 0], [0, 0, 0, 1.219, 0, 0], [0, 0, 0, 0, 1.219, 0], [0, 0, 0, 0, 0, 1.219]]
    G6 = [[0.0171364731454, 0, 0, 0, 0, 0], [0, 0.0171364731454, 0, 0, 0, 0], [0, 0, 0.033822, 0, 0, 0], [0, 0, 0, 0.1879, 0, 0], [0, 0, 0, 0, 0.1879, 0], [0, 0, 0, 0, 0, 0.1879]]

    Glist = [G1, G2, G3, G4, G5, G6]

    L1 = .425
    L2 = .39225
    W1 = .109
    H1 = .089159
    H2 = .09465
    Slist = [[0, 0, 1, 0, 0, 0], [0, 1, 0, -H1, 0, 0], [0, 1, 0, -H1, 0, L1], [0, 1, 0, -H1, 0, L1 + L2], [0, 0, -1, -W1, L1 + L2, 0], [0, 1, 0, H2 - H1, 0, L1 + L2]]

    M01 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, .089159,1]]).T
    M12 = np.array([[0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [.28, .13585, 0, 1]]).T
    M23 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, -.1197, .425-.03 , 1]]).T
    M34 = np.array([[0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0 , 0, .39225 - .25, 1]]).T
    M45 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, .093, 0, 1]]).T
    M56 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, .09465, 1]]).T

    Mlist = [M01, M12, M23, M34, M45, M56]

    Ftipmat = 0
    g = [0,0,-9.8]
    dt = 0.01
    intRes = 10

    simSec = 3
    taumat = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]] * int((simSec/dt))

    return ForwardDynamicsTrajectory(thetalist, dthetalist, taumat, g, Ftipmat, Mlist, Glist, Slist, dt, intRes)

def writeToFile(jointangles):
    f = open('fallingSim.csv','w')
    for i in range(0, len(jointangles)):
        joint = jointangles[i]
        f.write(",".join(map(str, joint)) + "\n")
    f.close()

def main():
    values = fallingSim()
    print "Complete!"
    thetamat = values[0]
    writeToFile(thetamat)
    theta1 = []
    theta2 = []
    theta3 = []
    theta4 = []
    theta5 = []
    theta6 = []
    for i in range(0, len(thetamat)):
        theta1.append(thetamat[i][0])
        theta2.append(thetamat[i][1])
        theta3.append(thetamat[i][2])
        theta4.append(thetamat[i][3])
        theta5.append(thetamat[i][4])
        theta6.append(thetamat[i][5])
    N = len(thetamat)
    Tf = len(thetamat) * .01
    timestamp = np.linspace(0, Tf, N)
    plt.plot(timestamp, theta1, label="Theta1")
    plt.plot(timestamp, theta2, label="Theta2")
    plt.plot(timestamp, theta3, label="Theta3")
    plt.plot(timestamp, theta4, label="Theta4")
    plt.plot(timestamp, theta5, label="Theta5")
    plt.plot(timestamp, theta6, label="Theta6")
    plt.ylim(-10, 10)
    plt.legend(loc='lower right')
    plt.xlabel("Time")
    plt.ylabel("Joint Angles/Velocities")
    plt.title("Plot of Joint Angles")
    plt.show()


if __name__=="__main__":
    main()
