module top(x0, x1, x2, x3, x4, y0);
input x0, x1, x2, x3, x4;
output y0;
wire t0, t1, t2, t3;
assign t0 = x0 & x4;
assign t1 = x1 & x4;
assign t2 = x2 & x4;
assign t3 = x3 & x4;
wire s1;
assign s1 = t0 ^ t1;
wire s2;
assign s2 = s1 ^ t2;
wire s3;
assign s3 = s2 ^ t3;
assign y0 = s3;
endmodule
