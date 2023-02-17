package main

import (
	"bullet/rangproof"
	"fmt"
	"math/big"
)

func main() {
	//ec := rangproof.GetECPrimeGroupKey()
	ran, com1, open1 := rangproof.RPProve(big.NewInt(5), nil)
	fmt.Println(ran)
	fmt.Println(com1)
	fmt.Println(open1)
	if rangproof.RPVerify(ran, com1) {
		fmt.Println("Range Proof Verification works")
	} else {
		fmt.Println("*****Range Proof FAILURE")
	}
	ranall, comall, _ := rangproof.SubProof(big.NewInt(5), big.NewInt(5), com1, open1)
	comv, _ := rangproof.PedersenSubNum(com1, big.NewInt(5))
	if comall.Comm.Equal(comv.Comm) {
		fmt.Println("commit right")
	}
	comvv := new(rangproof.PedersenCommit)
	comvv.Comm = comv.Comm
	if rangproof.RPVerify(ranall, comvv) {
		fmt.Println("Range Proof Verification works")
	} else {
		fmt.Println("*****Range Proof FAILURE")
	}

	//BigInteger number = BigInteger.valueOf(5);
	//BigInteger randomness = ProofUtils.randomNumber();
	//
	//GeneratorParams parameters = GeneratorParams.generateParams(256,curve);
	//GroupElement v = parameters.getBase().commit(number, randomness);
	//PeddersenCommitment<?> witness = new PeddersenCommitment<>(parameters.getBase(),number, randomness);
	//BouncyCastleECPoint.addCount=0;
	//BouncyCastleECPoint.expCount=0;
	//RangeProof proof = new RangeProofProver().generateProof(parameters, v, witness);
	//System.out.println(BouncyCastleECPoint.expCount);
	//System.out.println(BouncyCastleECPoint.addCount);
	//RangeProofVerifier verifier = new RangeProofVerifier();
	//verifier.verify(parameters, v, proof);

	//rp := new(RangeProof.BulletProof)
	//rp.Setup(4)
	//
	////生成commit
	////randomness, _ := rand.Int(rand.Reader, btcec.S256().N)
	////pd := new(RangeProof.PedersenBase)
	////V := pd.Commit(big.NewInt(5), randomness)
	//////witness
	////witness := new
	//
	//ran, com1, open1 := rp.Proof(big.NewInt(5), nil)
	//fmt.Println(ran)
	//fmt.Println(com1)
	//fmt.Println(open1)
	//
	//if rp.Verify(ran, com1) {
	//	fmt.Println("Range Proof Verification works")
	//} else {
	//	fmt.Println("*****Range Proof FAILURE")
	//}
	////com2, _ := rp.PedersenSubNum(com1, big.NewInt(5))
	////fmt.Println(com2)
	//ranall, comall, _ := rp.SubProof(big.NewInt(5), big.NewInt(5), com1, open1)
	//////_, comall, _ := rp.SubProof(big.NewInt(5), big.NewInt(5), com1, open1)
	////
	////给验证者com1,ran,5，ranall
	//comv, _ := rp.PedersenSubNum(com1, big.NewInt(5))
	//if comall.Comm.Equal(comv.Comm) {
	//	fmt.Println("commit right")
	//}
	//comvv := new(RangeProof.PedersenCommit)
	//comvv.Comm = comv.Comm
	//if rp.Verify(ranall, comvv) {
	//	fmt.Println("Range Proof Verification works")
	//} else {
	//	fmt.Println("*****Range Proof FAILURE")
	//}

	//把3个commit都上链，计算c* = c1-c2,如果c*==c，再验证c的proof

	//ran2 := rp.Proof(big.NewInt(0))
	//fmt.Println(ran2.Comm)

	//ran2.Comm = com2
	//if rp.Verify(ran2) {
	//	fmt.Println("Range Proof Verification works")
	//} else {
	//	fmt.Println("*****Range Proof FAILURE")
	//}
	//TODO:字符串转为结构体(编码方式，应该输出byte)
	//if rp.Verify(rp.Proof(big.NewInt(5))) {
	//	fmt.Println("Range Proof Verification works")
	//} else {
	//	fmt.Println("*****Range Proof FAILURE")
	//}

	//rangproof.EC = rangproof.NewECPrimeGroupKey(4) //输入n
	//// Testing smallest number in range
	//if rangproof.RPVerify(rangproof.RPProve(big.NewInt(3))) {
	//	fmt.Println("Range Proof Verification works")
	//} else {
	//	fmt.Println("*****Range Proof FAILURE")
	//}
}
