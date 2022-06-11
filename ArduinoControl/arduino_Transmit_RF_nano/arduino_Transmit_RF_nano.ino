/*
* Arduino Wireless Communication Tutorial
*     Example 1 - Transmitter Code
*                
* by Dejan Nedelkovski, www.HowToMechatronics.com
* 
* Library: TMRh20/RF24, https://github.com/tmrh20/RF24/
*/

#include <SPI.h>
#include <nRF24L01.h>
#include <RF24.h>

RF24 radio(7, 8); // CE, CSN

const byte address[6] = "00001";

typedef struct Packet
{
  uint32_t intensity;
  uint32_t vA;
  uint32_t vB;
  uint32_t vC;
  uint32_t vD;
  
  uint32_t PosX;
  uint32_t PosY;
};

Packet send_packet;


void setup() {
  Serial.begin(9600);
  radio.begin();
  radio.openWritingPipe(address);
  radio.setPALevel(RF24_PA_MIN);
  radio.stopListening();

  send_packet.intensity = 42;
  send_packet.PosX = 52;
  send_packet.PosY = 99;
}

void loop() {
  //const char text[] = "Hello World";
  
  //radio.write(&text, sizeof(text));
  radio.write(&send_packet, sizeof(Packet));
  delay(1000);
}
