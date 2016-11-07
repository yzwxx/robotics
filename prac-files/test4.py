import brickpi
import time

# Initializing Brickpi
interface=brickpi.Interface()
interface.initialize()



# Setting up sensor ports
port0 = 0 # port which ultrasoic sensor is plugged in to
port1 = 1
#port2 = 2

#interface.sensorEnable(port, brickpi.SensorType.SENSOR_ULTRASONIC);
#interface.sensorEnable(port, brickpi.SensorType.SENSOR_TOUCH);
#interface.sensorEnable(port0, brickpi.SensorType.SENSOR_TOUCH);
#interface.sensorEnable(port1, brickpi.SensorType.SENSOR_TOUCH);
#interface.sensorEnable(port2, brickpi.SensorType.SENSOR_ULTRASONIC);





motors = [0,1]

#print "input Kp:"
#k_p=input()
k_p = 480
k_i = 500 #1000
k_d = 5


def RotateAngle(angle):
    #angle = float(input("Enter a angle to rotate (in radians): "))
    interface.increaseMotorAngleReferences(motors,[-1.0*angle,angle])


    motorAngles = interface.getMotorAngles(motors)
    initialValues = [motorAngles[0][0],motorAngles[1][0]]
    
    while not interface.motorAngleReferencesReached(motors) :
        time.sleep(0.1)
    print "Destination reached!"
    

def GoLineAngle(angle):
    #angle = float(input("Enter a length to go (in radians): "))
    #interface.increaseMotorAngleReferences(motors,[1.006*angle,angle])
    
    #while True:
     #   usReading0 = interface.getSensorValue(port0)
      #  usReading1 = interface.getSensorValue(port1)   
    interface.increaseMotorAngleReferences(motors,[1.0*angle,angle])
    #print    (usReading0[0], usReading1[0])

    motorAngles = interface.getMotorAngles(motors)
    initialValues = [motorAngles[0][0],motorAngles[1][0]]
    
    
    
    while not interface.motorAngleReferencesReached(motors) :
        time.sleep(0.1)
    print "Destination reached!"



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
print motorParams0.pidParameters.K_d

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
a=0
b=0



#interface.startLogging("./log2_"+ str(k_p) +".txt")
LENGTH_M = 15   #15 FOR 40CM
LENGTH_N = LENGTH_M
# ANGLE = 16 17 19 18.50 
ANGLE = 20.5 #FOR 360

# Going Forward and Backwards

#GoLineAngle(LENGTH_M)
#time.sleep(0.1)
#GoLineAngle(-LENGTH_M)
#time.sleep(0.1)




# Turn Left and Turn Right
#RotateAngle(ANGLE)
#time.sleep(0.1)
#RotateAngle(-ANGLE)
#time.sleep(0.1)


# Drawing a square

#GoLineAngle(LENGTH_M)
time.sleep(0.1)
RotateAngle(ANGLE)
time.sleep(0.1)
#GoLineAngle(LENGTH_M)

'''
RotateAngle(ANGLE)
time.sleep(0.1)
GoLineAngle(LENGTH_N)
time.sleep(0.1)
RotateAngle(ANGLE)
time.sleep(0.1)
GoLineAngle(LENGTH_M)
time.sleep(0.1)
RotateAngle(ANGLE)
time.sleep(0.1)
GoLineAngle(LENGTH_N)
time.sleep(0.1)
RotateAngle(ANGLE)
time.sleep(0.1)
'''

#GoLineAngle(LENGTH_M) #for testing


#GoLineAngle(2*LENGTH_M)

#interface.stopLogging()
#fhandle.close
interface.terminate()





    


