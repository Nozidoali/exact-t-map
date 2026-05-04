module top(x0, x1, x2, y0);
input x0, x1, x2;
output y0;
wire n0;
assign n0 = x0 & x1;
assign y0 = n0 & x2;
endmodule
