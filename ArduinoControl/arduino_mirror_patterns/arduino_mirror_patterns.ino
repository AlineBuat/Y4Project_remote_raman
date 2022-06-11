#include "math.h"
#define RES 12


// pins used on the board
int QPpin_A = A0; //x = 1, y = 1
int QPpin_B = A1; // x = -1, y = 1;
int QPpin_C = A2; // x = -1, y = -1;
int QPpin_D = A3; // x = 1, y = -1

int mirror1Pin = DAC0;
int mirror2Pin = DAC1;

int blinkPin = LED_BUILTIN;
int blinkState = LOW;

unsigned long previousMicros = 0;
unsigned long currentMicros;


int maxNum = pow(2, RES)-1;
float freq = 5; // frequency in Hz
int T = 1000000/freq; // Period in microseconds


float currentPos1 = 0;
float currentPos2 = 0;

struct myDAC
{
  // relation between normal AnalogWrite() output and actual
  float p1 = 1; 
  int p2 = 0;
};

  myDAC myDAC0;
  myDAC myDAC1;

// enter calibration parameters for DAC0
float p1 = 1.1628;
float p2 = 532.4366;


void setup() {
  // put your setup code here, to run once:
  pinMode(blinkPin, OUTPUT);
  pinMode(mirror1Pin, OUTPUT);
  analogWriteResolution(RES);
  pinMode(mirror2Pin, OUTPUT);

  analogWriteResolution(RES);
  analogReadResolution(RES);

  SerialUSB.begin(9600);
  while(!SerialUSB){}
  SerialUSB.print("MaxNum: ");
  SerialUSB.println(maxNum);

  // initialize DACs offsets
  myDAC0.p1 = 1.1628;
  myDAC0.p2 = 532.4366;
  myDAC1.p1 = 0.9649;
  myDAC1.p2 = 68.25;
}

void loop() {
  // put your main code here, to run repeatedly:

  float weight;
  float blinkSignal;
  float mirror1Sig;
  float mirror2Sig;

  currentMicros = micros();
  
  while (SerialUSB.available() > 0)
  {
    float posWanted1 = SerialUSB.parseFloat();
    float posWanted2 = SerialUSB.parseFloat(); // read second valid float
    
    currentPos1 = posWanted1;
    currentPos2 = posWanted2;
    SerialUSB.print("I received a bytyMcByte Face!  ");
    SerialUSB.println(currentPos1, DEC);
    SerialUSB.println(currentPos2, DEC);
    mirror1Sig = (currentPos1 + 1)*maxNum/2;
    mirror2Sig = (currentPos2 + 1)*maxNum/2;

    SerialUSB.print("New mirror signal: ");
    SerialUSB.print("Mirror1: ");
    SerialUSB.print(mirror1Sig);
    SerialUSB.print(" , Mirror2: ");
    SerialUSB.println(mirror2Sig);
    
  }
   

  
  mirror2Sig = f_test_triangle(currentMicros, T, maxNum);
  mirror1Sig = f_test_triangle(currentMicros, T, maxNum);

//  mirror1Sig = (currentPos1 + 1)* maxNum / 2;
//  mirror2Sig = (currentPos2 + 1)* maxNum / 2;
  
  // correct the output range
  mirror1Sig = myDAC0.p1*mirror1Sig + myDAC0.p2;
  mirror2Sig = myDAC1.p1*mirror2Sig + myDAC1.p2;

  // correct for saturating values (not higher than maximal analog authorised values)
  if (mirror1Sig > maxNum) mirror1Sig = maxNum;
  else if (mirror1Sig <0) mirror1Sig = 0;
  if (mirror2Sig > maxNum) mirror2Sig = maxNum;
  else if (mirror2Sig <0) mirror2Sig = 0;
  
  analogWrite(mirror1Pin, mirror1Sig);
  analogWrite(mirror2Pin, mirror2Sig);
  
  //analogWrite(blinkPin, blinkSignal);

  if (currentMicros - previousMicros > 1000000)
  {
    blinkState = !blinkState;
    previousMicros = currentMicros;
    SerialUSB.println("Here are the signals!");
    SerialUSB.println(mirror1Sig);
    SerialUSB.println(mirror2Sig);
    digitalWrite(blinkPin, blinkState);
  }
  
}


float f_test_sin(float t, float T, int max_val)
{ 
  /***function to return. Here a sinusoid of period T, amplitude N, centered around N/2. Returns an int***/
  int val = max_val/2*sin(2*PI/T*t) + max_val/2;
  
  return val;
}

float f_test_triangle(float t, float T, int max_val)
{
  /*Function returns a triangle function between 0 and max_val, and period T*/
  int val = 0;
  float a = 2*max_val/T; // how fast does our ramp grow
  t = ((int) t) % ((int)T); // remainder of time by period T
  if ( t < T/2)
  {
    val = a*t;
  }
  else // case we are between T/2 and T
  {
    val = a*(T - t);
  }
  return val;
}

float f_test_ramp(float t, float T, int max_val)
{
  int val = 0;
  float a = max_val/T;
  t = ((int) t) % ((int)T); // remainder of time by period T
  
  val = t*a;
  return val;
}

float f_test_square(float t, float T, int max_val)
{
  int val = 0;
  t = ((int) t) % ((int)T); // remainder of time by period T
  if (t> T/2)
  {
    val = max_val;
  } else
  {
    val = 0;
  }
  return val;
}
