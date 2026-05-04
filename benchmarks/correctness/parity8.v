module top(x0, x1, x2, x3, x4, x5, x6, x7, y0);
input x0, x1, x2, x3, x4, x5, x6, x7;
output y0;
wire t0, t1, t2, t3, t4, t5;
assign t0 = x0 ^ x1;
assign t1 = x2 ^ x3;
assign t2 = x4 ^ x5;
assign t3 = x6 ^ x7;
assign t4 = t0 ^ t1;
assign t5 = t2 ^ t3;
assign y0 = t4 ^ t5;
endmodule
