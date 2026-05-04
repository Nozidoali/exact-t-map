module top(a0, a1, b0, b1, cin, s0, s1, cout);
input a0, a1, b0, b1, cin;
output s0, s1, cout;
wire ax0, ac0, t0, c1, ax1, ac1, t1;
assign ax0 = a0 ^ b0;
assign ac0 = a0 & b0;
assign s0 = ax0 ^ cin;
assign t0 = ax0 & cin;
assign c1 = ac0 | t0;
assign ax1 = a1 ^ b1;
assign ac1 = a1 & b1;
assign s1 = ax1 ^ c1;
assign t1 = ax1 & c1;
assign cout = ac1 | t1;
endmodule
