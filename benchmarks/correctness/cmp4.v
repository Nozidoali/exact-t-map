module top(a0, a1, a2, a3, b0, b1, b2, b3, y0);
input a0, a1, a2, a3, b0, b1, b2, b3;
output y0;
wire d0, d1, d2, d3;
wire e0, e1, e2, e3;
wire p0, p1;
assign d0 = a0 ^ b0;
assign d1 = a1 ^ b1;
assign d2 = a2 ^ b2;
assign d3 = a3 ^ b3;
assign e0 = ~d0;
assign e1 = ~d1;
assign e2 = ~d2;
assign e3 = ~d3;
assign p0 = e0 & e1;
assign p1 = e2 & e3;
assign y0 = p0 & p1;
endmodule
