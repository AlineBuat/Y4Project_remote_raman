/*
* Arduino Wireless Communication Tutorial
*       Example 1 - Receiver Code
*                
* by Dejan Nedelkovski, www.HowToMechatronics.com
* 
* Library: TMRh20/RF24, https://github.com/tmrh20/RF24/
*/

#include <SPI.h>
#include <nRF24L01.h>
#include <RF24.h>

#define RES_RECEIVE 10;

RF24 radio(7, 8); // CE, CSN

const byte address[6] = "00001";

struct Packet
{
  uint32_t intensity;
  uint32_t vA;
  uint32_t vB;
  uint32_t vC;
  uint32_t vD;
  
  int32_t PosX;
  int32_t PosY;
};


unsigned long previous_millis = 0;
int blinkState = LOW;

Packet read_packet;

void setup() {
  pinMode(LED_BUILTIN, OUTPUT);
  
  SerialUSB.begin(9600);
  while (!SerialUSB){}
  SerialUSB.println("Program Started!");
  radio.begin();
  radio.openReadingPipe(0, address);
  radio.setPALevel(RF24_PA_MIN);

  radio.startListening();
}
void loop() {
  unsigned long current_millis = millis();
  
  if (radio.available()) {
//    SerialUSB.println("Message Received!");
    //char text[32] = "";
    //radio.read(&text, sizeof(text));
    radio.read(&read_packet, sizeof(Packet));
    
//    SerialUSB.print(read_packet.PosX);
//    SerialUSB.print('\t');
//    SerialUSB.print(read_packet.PosY);
    SerialUSB.print("vA");
    SerialUSB.print(read_packet.vA);
    SerialUSB.print('\t');
    SerialUSB.print("vB");
    SerialUSB.print(read_packet.vB);
    SerialUSB.print('\t');
    SerialUSB.print("vC");
    SerialUSB.print(read_packet.vC);
    SerialUSB.print('\t');
    SerialUSB.print("vD");
    SerialUSB.print(read_packet.vD);
    SerialUSB.print('\n');
  }

  if (current_millis - previous_millis >2000)
  {
    previous_millis = current_millis;
    blinkState = !blinkState;
    digitalWrite(LED_BUILTIN, blinkState);
//    SerialUSB.println("I'm still running!");
  }

  delay(1);
}
