import brickpi
import time
import math

robotX = 0
robotY = 0
robotTheta = 0

interface = brickpi.Interface()
interface.initialize()
motors = [0, 1]

k_p = 480
k_i = 400
k_d = 5

ANGLEPRECETIMETER = 11.9493 / 40
ANGLEPERDEGREE = 4.15 / 90

left_coefficient = 1.035

def rotate(rotation):
    angle = rotation * ANGLEPERDEGREE
    interface.increaseMotorAngleReferences(motors, [-left_coefficient * angle, angle])
    motorAngles = interface.getMotorAngles(motors)
    initialValues = [motorAngles[0][0], motorAngles[1][0]]
    while not interface.motorAngleReferencesReached(motors):
        time.sleep(0.1)
    time.sleep(0.1)

def goLine(distance):
    angle = distance * ANGLEPRECETIMETER
    interface.increaseMotorAngleReferences(motors, [left_coefficient * angle, angle])
    motorAngles = interface.getMotorAngles(motors)
    initialValues = [motorAngles[0][0], motorAngles[1][0]]
    while not interface.motorAngleReferencesReached(motors):
        time.sleep(0.1)
    time.sleep(0.1)

def navigateToWaypoint(x, y):
	global robotX
	global robotY
	global robotTheta
	dx = x - robotX
	dy = y - robotY
	beta = math.atan2(dy, dx) / math.pi * 180
	dtheta = robotTheta - beta
	if dtheta > 180:
		dtheta = dtheta - 360
	rotate(dtheta)
	robotTheta = beta
	distance = math.sqrt(dx * dx + dy * dy)
	goLine(distance)
	robotX = x
	robotY = y

#Initialization
interface.motorEnable(motors[0])
interface.motorEnable(motors[1])

motorParams0 = interface.MotorAngleControllerParameters()
motorParams0.maxRotationAcceleration = 6
motorParams0.maxRotationSpeed = 12
motorParams0.feedForwardGain = 255/20.0
motorParams0.minPWM = 18.0
motorParams0.pidParameters.minOutput = -255
motorParams0.pidParameters.maxOutput = 255
motorParams0.pidParameters.k_p = k_p
motorParams0.pidParameters.k_i = k_i
motorParams0.pidParameters.K_d = k_d

motorParams1 = interface.MotorAngleControllerParameters()
motorParams1.maxRotationAcceleration = 6.0
motorParams1.maxRotationSpeed = 12
motorParams1.feedForwardGain = 255/20.0
motorParams1.minPWM = 18.0
motorParams1.pidParameters.minOutput = -255
motorParams1.pidParameters.maxOutput = 255
motorParams1.pidParameters.k_p = k_p
motorParams1.pidParameters.k_i = k_i
motorParams1.pidParameters.K_d = k_d


interface.setMotorAngleControllerParameters(motors[0],motorParams0)
interface.setMotorAngleControllerParameters(motors[1],motorParams1)

#interface.startLogging("./log2_" + str(k_p) + ".txt"

navigateToWaypoint(-10,0)

#interface.stopLogging()
interface.terminate()
