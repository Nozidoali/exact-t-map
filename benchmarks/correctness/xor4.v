module top(x0, x1, x2, x3, y0);
input x0, x1, x2, x3;
output y0;
wire n0, n1;
assign n0 = x0 ^ x1;
assign n1 = x2 ^ x3;
assign y0 = n0 ^ n1;
endmodule
