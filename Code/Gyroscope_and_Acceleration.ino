
#include <Wire.h>
#include <SD.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>
#include "RTClib.h"
Adafruit_BNO055 bno = Adafruit_BNO055();
RTC_DS1307 RTC;

#define SAMPLE_PERIOD 10
#define WRITE_PERIOD 10*SAMPLE_PERIOD
#define SAVE_MAGNETOMETER_DATA 0
#define SAVE_GYROSCOPE_DATA 1
#define SAVE_EULER_ANGLE_DATA 1
#define SAVE_ACCELEROMETER_DATA 1
#define SAVE_LIN_ACCELERATION_DATA 1
#define SAVE_GRAVITY_DATA 0
#define ECHO_TO_SERIAL 1
#define WAIT_TO_START 0

#define redLEDpin 3
#define greenLEDpin 2

uint32_t writeTime = 0;
const int chipSelect = 10;

File logfile;

void setup(void) {
  // put your setup code here, to run once:
  Serial.begin(9600);
  Serial.println();
  initialize_leds();
  initialize_sd_card();
  Wire.begin();
  initialize_rtc();
  
  initialize_bno055();  
  create_data_header();
  
}

void loop(void) {
  // put your main code here, to run repeatedly:
  DateTime now;
  delay((SAMPLE_PERIOD - 1) - (millis() % SAMPLE_PERIOD));
  read_and_write_sensor_data();

  if((millis() - writeTime) < WRITE_PERIOD){
    return;
  }
  writeTime = millis();
  logfile.flush();
  
  
}


void initialize_sd_card(void){
  Serial.print("Initializing SD card --- ");
  pinMode(chipSelect, OUTPUT);

  if(!SD.begin(chipSelect)){
    error("SD card failed or not present");
  }

  Serial.println("SD card initialized");

  char filename[] = "LOGGER00.CSV";
  for (uint8_t i = 0; i< 100; i++){
    filename[6] = i/10 + '0';
    filename[7] = i%10 + '0';
    if(!SD.exists(filename)){  
      logfile = SD.open(filename, FILE_WRITE);
      break;
    }
  }

  if(!logfile){
    error("Could not create data file on SD card."); 
  }
 
  Serial.print("Logging to: ");
  Serial.println(filename);
}

void initialize_rtc(void){
  if(!RTC.begin()){
    error("RTC communication failed!");
  }
}

void initialize_leds(void){
  pinMode(redLEDpin,OUTPUT);
  pinMode(greenLEDpin, OUTPUT);
  digitalWrite(redLEDpin, LOW);
  digitalWrite(greenLEDpin, LOW);
}

void initialize_bno055(void){
  print_SD_and_serial("BNO Section");
  logfile.print("BNO Section");
  if(!bno.begin()){
    error("Cannot connect to BNO055 board");
  }
  bno.setExtCrystalUse(true);
}

void error(char const *str){
  Serial.print("error: ");
  Serial.println(str);
  while(1){
    digitalWrite(redLEDpin, HIGH);
    digitalWrite(greenLEDpin, HIGH);
    delay(500);
    digitalWrite(redLEDpin, LOW);
    digitalWrite(greenLEDpin, LOW);
    delay(500);
  }
}

void print_SD_and_serial(String str){
  logfile.print(str);
  #if ECHO_TO_SERIAL
    Serial.print(str);
  #endif
}
void println_SD_and_serial(String str){
  logfile.println(str);
  #if ECHO_TO_SERIAL
    Serial.print(str);
  #endif
}

void display_calibration_status(uint8_t cal_system, 
                                uint8_t cal_magnetometer, 
                                uint8_t cal_gyroscope, 
                                uint8_t cal_accelerometer){ 
  uint8_t cal_minimum = 3;
  if (cal_system < cal_minimum){
    cal_minimum = cal_system;
  }

  #if (SAVE_MAGNETOMETER_DATA)
    if(cal_magnetometer <cal_minimum){
      cal_minimum = cal_magnetometer;
    }
  #endif

  #if ((SAVE_GYROSCOPE_DATA) || (SAVE_EULER_ANGLE_DATA))
    if(cal_gyroscope <cal_minimum){
      cal_minimum = cal_gyroscope;
    }
  #endif
  #if ((SAVE_ACCELEROMETER_DATA)||(SAVE_LIN_ACCELERATION_DATA)||(SAVE_GRAVITY_DATA))
    if(cal_accelerometer <cal_minimum){
      cal_minimum = cal_accelerometer;
    }
  #endif

  
  if (cal_minimum == 0){
    digitalWrite(redLEDpin, LOW);
    digitalWrite(greenLEDpin, LOW);
  }
  if (cal_minimum == 1){
    digitalWrite(redLEDpin, HIGH);
    digitalWrite(greenLEDpin, LOW);
  }
  if (cal_minimum == 2){
    digitalWrite(redLEDpin, LOW);
    digitalWrite(greenLEDpin, HIGH);
  }

  if (cal_minimum == 3){
    digitalWrite(redLEDpin, HIGH);
    digitalWrite(greenLEDpin, HIGH);
  }
}

void create_data_header(void){
  DateTime rtc_time_1 = RTC.now();
  print_SD_and_serial(String("Unix Time: "));
  print_SD_and_serial(String(rtc_time_1.unixtime()));
  println_SD_and_serial(String(""));
  print_SD_and_serial(String("Standard Time: "));
    print_SD_and_serial(String('"'));
  print_SD_and_serial(String(rtc_time_1.year()));
  print_SD_and_serial(String("/"));
  print_SD_and_serial(String(rtc_time_1.month()));
  print_SD_and_serial(String("/"));
  print_SD_and_serial(String(rtc_time_1.day()));
  print_SD_and_serial(String(" "));
  print_SD_and_serial(String(rtc_time_1.hour()));
  print_SD_and_serial(String(":"));
  print_SD_and_serial(String(rtc_time_1.minute()));
  print_SD_and_serial(String(":"));
  print_SD_and_serial(String(rtc_time_1.second()));
  print_SD_and_serial(String('"'));
  println_SD_and_serial(String(""));


  print_SD_and_serial(String("ms"));
  #if (SAVE_MAGNETOMETER_DATA)
    print_SD_and_serial(String(", uT, uT, uT"));
  #endif
  #if (SAVE_GYROSCOPE_DATA)
    print_SD_and_serial(String(", 1/s, 1/s, 1/s"));
  #endif
  #if (SAVE_EULER_ANGLE_DATA)
    print_SD_and_serial(String(", degree, degree, degree"));
  #endif
  #if (SAVE_ACCELEROMETER_DATA)
    print_SD_and_serial(String(", m/s^2, m/s^2, m/s^2"));
  #endif
  #if (SAVE_LIN_ACCELERATION_DATA)
    print_SD_and_serial(String(", m/s^2, m/s^2, m/s^2"));
  #endif
  #if (SAVE_GRAVITY_DATA)
    print_SD_and_serial(String(", m/s^2, m/s^2, m/s^2"));
  #endif
  print_SD_and_serial(String(", N/A"));
  #if (SAVE_MAGNETOMETER_DATA)
    print_SD_and_serial(String(", N/A"));
  #endif
  #if (SAVE_GYROSCOPE_DATA||SAVE_EULER_ANGLE_DATA)
    print_SD_and_serial(String(", N/A"));
  #endif
  #if (SAVE_ACCELEROMETER_DATA||SAVE_LIN_ACCELERATION_DATA||SAVE_GRAVITY_DATA)
   print_SD_and_serial(String(", N/A"));
  #endif
  println_SD_and_serial(String(""));

  print_SD_and_serial(String("t"));
  #if (SAVE_MAGNETOMETER_DATA)
   print_SD_and_serial(String(", B_x, B_y, B_z"));
  #endif
  #if (SAVE_GYROSCOPE_DATA)    
    print_SD_and_serial(String(", omega_x, omega_y, omega_z"));
  #endif
  #if (SAVE_EULER_ANGLE_DATA) 
    print_SD_and_serial(String(", theta_x, theta_y, theta_z"));
  #endif
  #if (SAVE_ACCELEROMETER_DATA)
   print_SD_and_serial(String(", a_x, a_y, a_z"));
  #endif
  #if (SAVE_LIN_ACCELERATION_DATA) 
    print_SD_and_serial(String(", a_lin_x, a_lin_y, a_lin_z"));
  #endif
  #if (SAVE_GRAVITY_DATA)
    print_SD_and_serial(String(", a_g_x, a_g_y, a_g_z"));
  #endif
  print_SD_and_serial(String(", cal_system"));
  #if (SAVE_MAGNETOMETER_DATA)
   print_SD_and_serial(String(", cal_magnet"));
  #endif
  #if (SAVE_GYROSCOPE_DATA||SAVE_EULER_ANGLE_DATA)
    print_SD_and_serial(String(", cal_gyro"));
  #endif
  #if (SAVE_ACCELEROMETER_DATA||SAVE_LIN_ACCELERATION_DATA||SAVE_GRAVITY_DATA)
   print_SD_and_serial(String(", cal_accel"));
  #endif
  println_SD_and_serial(String(""));
   
}

void read_and_write_sensor_data(void){
  #if (SAVE_MAGNETOMETER_DATA)
    imu::Vector<3> magnetometer_vector;
  #endif
  #if (SAVE_GYROSCOPE_DATA)
    imu::Vector<3> gyroscope_vector;
  #endif
  #if (SAVE_EULER_ANGLE_DATA)
    imu::Vector<3> euler_vector;
  #endif
  #if (SAVE_ACCELEROMETER_DATA)
    imu::Vector<3> accelerometer_vector;
  #endif
  #if (SAVE_LIN_ACCELERATION_DATA)
    imu::Vector<3> linearaccel_vector;
  #endif
  #if (SAVE_GRAVITY_DATA)
    imu::Vector<3> gravity_vector;
  #endif
  
  uint8_t cal_system, cal_magnetometer, cal_gyroscope, cal_accelerometer;
  cal_system = 0;
  cal_magnetometer = 0;
  cal_gyroscope = 0;
  cal_accelerometer = 0;

  bno.getCalibration(&cal_system, &cal_gyroscope, &cal_accelerometer, &cal_magnetometer);
  display_calibration_status( cal_system, cal_magnetometer, cal_gyroscope, cal_accelerometer );

  #if (SAVE_MAGNETOMETER_DATA)
    magnetometer_vector = bno.getVector(Adafruit_BNO055::VECTOR_MAGNETOMETER);
  #endif
  #if (SAVE_GYROSCOPE_DATA)
    gyroscope_vector = bno.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);
  #endif
  #if (SAVE_EULER_ANGLE_DATA)
    euler_vector = bno.getVector(Adafruit_BNO055::VECTOR_EULER);
  #endif
  #if (SAVE_ACCELEROMETER_DATA)
    accelerometer_vector = bno.getVector(Adafruit_BNO055::VECTOR_ACCELEROMETER);
  #endif  
  #if (SAVE_LIN_ACCELERATION_DATA)
    linearaccel_vector = bno.getVector(Adafruit_BNO055::VECTOR_LINEARACCEL);
  #endif
  #if (SAVE_GRAVITY_DATA)
    gravity_vector = bno.getVector(Adafruit_BNO055::VECTOR_GRAVITY);
  #endif
  
  unsigned long mcr_time_ms = millis();
  print_SD_and_serial(String(mcr_time_ms));
  #if (SAVE_MAGNETOMETER_DATA)
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(magnetometer_vector.x()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(magnetometer_vector.y()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(magnetometer_vector.z()));
  #endif
  #if (SAVE_GYROSCOPE_DATA)
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(gyroscope_vector.x()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(gyroscope_vector.y()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(gyroscope_vector.z()));
  #endif
  #if (SAVE_EULER_ANGLE_DATA) 
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(euler_vector.x()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(euler_vector.y()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(euler_vector.z()));
  #endif
  #if (SAVE_ACCELEROMETER_DATA)
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(accelerometer_vector.x()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(accelerometer_vector.y()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(accelerometer_vector.z()));
  #endif
  #if (SAVE_LIN_ACCELERATION_DATA) 
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(linearaccel_vector.x()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(linearaccel_vector.y()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(linearaccel_vector.z()));
  #endif
  #if (SAVE_GRAVITY_DATA)
    // This writes the gravitational acceleration vector component variables
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(gravity_vector.x()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(gravity_vector.y()));
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(gravity_vector.z()));
  #endif
  print_SD_and_serial(String(", "));
  print_SD_and_serial(String(cal_system));
  #if (SAVE_MAGNETOMETER_DATA)
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(cal_magnetometer));
  #endif
  #if (SAVE_GYROSCOPE_DATA||SAVE_EULER_ANGLE_DATA)
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(cal_gyroscope));
  #endif
  #if (SAVE_ACCELEROMETER_DATA||SAVE_LIN_ACCELERATION_DATA||SAVE_GRAVITY_DATA)
    print_SD_and_serial(String(", "));
    print_SD_and_serial(String(cal_accelerometer));
  #endif

  println_SD_and_serial(String(""));
  
}

