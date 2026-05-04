module top(a0, a1, b0, b1, p0, p1, p2, p3);
input a0, a1, b0, b1;
output p0, p1, p2, p3;
wire pp00, pp01, pp10, pp11;
wire c1;
assign pp00 = a0 & b0;
assign pp01 = a0 & b1;
assign pp10 = a1 & b0;
assign pp11 = a1 & b1;
assign p0 = pp00;
assign p1 = pp01 ^ pp10;
assign c1 = pp01 & pp10;
assign p2 = pp11 ^ c1;
assign p3 = pp11 & c1;
endmodule
