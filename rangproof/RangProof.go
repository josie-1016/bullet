package rangproof

//
//import (
//	"math/big"
//)
//
//type BulletProof struct {
//	ECP *CryptoParams
//	rp  *RangeProof
//}
//
//func (bp *BulletProof) Setup(n int) {
//	bp.rp = new(RangeProof)
//	bp.ECP = new(CryptoParams)
//	bp.ECP.NewECPrimeGroupKey(n)
//}
//
//func (bp *BulletProof) Proof(v *big.Int, open *big.Int) (*RangeProof, *PedersenCommit, *big.Int) {
//	return bp.rp.RPProve(v, open, bp.ECP)
//}
//
//func (bp *BulletProof) Verify(rproof *RangeProof, commit *PedersenCommit) bool {
//	return rproof.RPVerify(commit, bp.ECP)
//	//return RPVerify(rproof, commit)
//}
//
//func (bp *BulletProof) SubProof(v1 *big.Int, v2 *big.Int, com1 *PedersenCommit, open1 *big.Int) (*RangeProof, *PedersenCommit, *big.Int) {
//	v := new(big.Int)
//	//commit,_ := bp.PedersenAddNum(com1.Comm, v2)
//	proof, com, open := bp.Proof(v.Sub(v1, v2), open1)
//	return proof, com, open
//}
//
//func (bp *BulletProof) PedersenAddNum(commit *PedersenCommit, v *big.Int) (*PedersenCommit, error) {
//	//com, err := bp.CommitTo(v)
//
//	//if err != nil {
//	//	return ECPoint{}, err
//	//}
//	//致盲因子不变
//	return commit.PedersenAddNum(v, bp.ECP)
//	//return bp.PedersenAddCommitment(cX, com), nil
//}
//
//func (bp *BulletProof) PedersenSubNum(commit *PedersenCommit, v *big.Int) (*PedersenCommit, error) {
//	//com, err := bp.CommitTo(v)
//	//if err != nil {
//	//	return ECPoint{}, err
//	//}
//	//return bp.PedersenSubCommitment(cX, com), nil
//	return commit.PedersenSubNum(v, bp.ECP)
//	//return cX.Add(EC.G.Mult(v).Neg()), nil
//}
