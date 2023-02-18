package main

import (
	"bullet/rangproof"
	"fmt"
	"math/big"
)

func main() {
	// 生成范围证明:eg.x>5
	//rangproof.EC = rangproof.NewECPrimeGroupKey(3)
	ran, com1, open1 := rangproof.RPProve(big.NewInt(5), nil) //参数为私密值和选择的致盲因子（可为空），返回证明、pedersen承诺和使用的致盲因子
	fmt.Println(ran)
	fmt.Println(com1)
	fmt.Println(open1)

	//验证证明，输入证明和pedersen承诺
	if rangproof.RPVerify(ran, com1) {
		fmt.Println("Range Proof Verification works")
	} else {
		fmt.Println("Range Proof FAILURE")
	}

	//参数为v1的承诺和使用的致盲因子，生成v1 - v2为非负整数的证明、v1 - v2的承诺和致盲因子
	ran2, com2, _ := rangproof.SubProof(big.NewInt(5), big.NewInt(5), com1, open1)

	//可以将生成的com1、com2和ran2发送给另一方，证明自己的值大于5
	//另一方先判断承诺是否正确，然后验证两者差值的证明，如通过则证明大于0，满足条件

	//可以进行“承诺 - 数值”的同态运算，致盲因子不变
	comv, _ := rangproof.PedersenSubNum(com1, big.NewInt(5))
	//测试是否与SubProof方法计算出的承诺值相同
	if com2.Comm.Equal(comv.Comm) {
		fmt.Println("commit right")
	}
	//验证
	if rangproof.RPVerify(ran2, comv) {
		fmt.Println("Range Proof Verification works")
	} else {
		fmt.Println("*****Range Proof FAILURE")
	}

}
