import brickpi
import time

interface=brickpi.Interface()
interface.initialize()

port0 = 0 # port which ultrasoic sensor is plugged in to
port1 = 1
port2 = 2
#interface.sensorEnable(port, brickpi.SensorType.SENSOR_ULTRASONIC);
#interface.sensorEnable(port, brickpi.SensorType.SENSOR_TOUCH);
interface.sensorEnable(port0, brickpi.SensorType.SENSOR_TOUCH);
interface.sensorEnable(port1, brickpi.SensorType.SENSOR_TOUCH);
interface.sensorEnable(port2, brickpi.SensorType.SENSOR_ULTRASONIC);


while True:
    usReading0 = interface.getSensorValue(port0)
    usReading1 = interface.getSensorValue(port1)
    usReading2 = interface.getSensorValue(port2)
    
    #if (usReading0 and usReading1 and usReading2):
       # print usReading0, usReading1, usReading2
    print usReading0,usReading1,usReading2
    
    
    #else:
     #   print "Failed US reading"

    time.sleep(0.2)

interface.terminate()