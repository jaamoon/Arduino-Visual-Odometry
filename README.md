# Arduino-Visual-Odometry

## Input

It receives a set of eight corresponding points via serial communication.
With the format u1,v1,u2,v2,...,u8,v8 \n.
After receiving the second set of eight points, it calculates the essential matrix and then the rotation matrix and the translation vector.

## Output

Using interrupts guarantees that every hundred milliseconds it outputs the comma-separated rotation matrix and the translation vector through the serial port.
