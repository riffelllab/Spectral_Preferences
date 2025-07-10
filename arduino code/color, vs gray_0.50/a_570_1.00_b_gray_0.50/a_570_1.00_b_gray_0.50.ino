#include <Wire.h>
#include <Adafruit_PWMServoDriver.h>

//defines pwm brightness (0-4096 with 0 off and 4096 fully on)
//for LEDs between 385-740 nm for LED box A and LED box B

// LED box A, pwm values for LEDs 385 - 740 nm
int a_385 = 0;
int a_395 = 0;
int a_415 = 0;
int a_430 = 0;
int a_450 = 0;
int a_470 = 0;
int a_505 = 0;
int a_525 = 0;
int a_545 = 0;
int a_570 = 2990;
int a_590 = 0;
int a_625 = 0;
int a_645 = 0;
int a_660 = 0;
int a_680 = 0;
int a_700 = 0;
int a_740 = 0;

// LED box B, pwm values for LEDs 385 - 740 nm
int b_385 = 0;
int b_395 = 0;
int b_415 = 0;
int b_430 = 0;
int b_450 = 828;
int b_470 = 0;
int b_505 = 440;
int b_525 = 743;
int b_545 = 949;
int b_570 = 608;
int b_590 = 86;
int b_625 = 875;
int b_645 = 0;
int b_660 = 0;
int b_680 = 0;
int b_700 = 0;
int b_740 = 0;

// defines an array of pwm values for LED box B, board 0 (pwm[0], 0x40)
// 0s in the list indicate unsed pins
int board_0[16] = { a_385,  a_470,  a_545,  a_570,
                    a_645,  a_660,  a_700,  0,
                    0,      0,      0,      a_395,
                    a_450,  a_570,  a_590,  a_625
                  };

// defines an array of pwm values for LED box B, board 1 (pwm[1], 0x41)
// 0s in the list indicate unsed pins
int board_1[16] = { a_385,  a_430,  a_545,  a_570, 
                    a_645,  a_680,  a_700,  0,
                    0,      0,      a_415,  a_505,
                    a_525,  a_570,  a_590,  a_740
                  };

// defines an array of pwm values for LED box B, board 2 (pwm[2], 0x43)
// 0s in the list indicate unsed pins
int board_2[16] = { b_385,  b_470,  b_545,  b_570,
                    b_645,  b_660,  b_700,  0,
                    0,      0,      0,      b_395,
                    b_450,  b_570,  b_590,  b_625
                  };

// defines an array of pwm values for LED box B, board 3 (pwm[3], 0x42)
// 0s in the list indicate unsed pins
int board_3[16] = { b_385,  b_430,  b_545,  b_570, 
                    b_645,  b_680,  b_700,  0,
                    0,      0,      b_415,  b_505,
                    b_525,  b_570,  b_590,  b_740
                  };

// defines the addresses of the four PCA9685 boards   
Adafruit_PWMServoDriver pwm[] =
  { // boards for LED box A
    Adafruit_PWMServoDriver(0x40),
    Adafruit_PWMServoDriver(0x41),
    // boards for LED box B
    Adafruit_PWMServoDriver(0x43),
    Adafruit_PWMServoDriver(0x42)
  };

void setup() {
  for (int b = 0; b <= 3; b++) {
    // sets up the I2C interface and hardware
    pwm[b].begin();
    // sets the pwm frequency
    pwm[b].setPWMFreq(1600);
    // sets the pins to open drain, pwm pin will be sink and V+ source
    pwm[b].setOutputMode(false);
  }

  //sets pwm pin values for LED box A, pwm[0]
  for (int i = 0; i <= 15; i++) {
    pwm[0].setPin(i,board_0[i],true);
  }

  //sets pwm pin values for LED box A, pwm[1]
  for (int i = 0; i <= 15; i++) {
    pwm[1].setPin(i,board_1[i],true);
  }


  //sets pwm pin values for LED box B, pwm[2]
  for (int i = 0; i <= 15; i++) {
    pwm[2].setPin(i,board_2[i],true);
  }

  //sets pwm pin values for LED box B, pwm[3]
  for (int i = 0; i <= 15; i++) {
    pwm[3].setPin(i,board_3[i],true);
  }
  
}

void loop() {
}
