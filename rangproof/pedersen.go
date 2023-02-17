package rangproof

import "math/big"

type PedersenCommit struct {
	//G ECPoint // G value for commitments of a single value
	//H ECPoint // H value for commitments of a single value
	Comm ECPoint
}

//func (p *PedersenBase) Commit(v *big.Int, r *big.Int) ECPoint {
//	return EC.G.Mult(v).Add(EC.H.Mult(r))
//}

//func (p *PedersenCommit) PedersenSubNum(v *big.Int, ec *CryptoParams) (*PedersenCommit, error) {
//	commit := new(PedersenCommit)
//	commit.Comm = p.Comm.Add(ec.G.Mult(v, ec).Neg(ec), ec)
//	return commit, nil
//}
//
//func (p *PedersenCommit) PedersenAddNum(v *big.Int, ec *CryptoParams) (*PedersenCommit, error) {
//	commit := new(PedersenCommit)
//	commit.Comm = p.Comm.Add(ec.G.Mult(v, ec), ec)
//	return commit, nil
//}
//
//func (p *PedersenCommit) PedersenAddCommit(com *PedersenCommit, ec *CryptoParams) (*PedersenCommit, error) {
//	commit := new(PedersenCommit)
//	commit.Comm = p.Comm.Add(com.Comm, ec)
//	return commit, nil
//}
//
//func (p *PedersenCommit) PedersenSubCommit(com *PedersenCommit, ec *CryptoParams) (*PedersenCommit, error) {
//	commit := new(PedersenCommit)
//	commit.Comm = p.Comm.Add(com.Comm.Neg(ec), ec)
//	return commit, nil
//}

func PedersenSubNum(p *PedersenCommit, v *big.Int) (*PedersenCommit, error) {
	commit := new(PedersenCommit)
	commit.Comm = p.Comm.Add(EC.G.Mult(v).Neg())
	return commit, nil
}

func PedersenAddNum(p *PedersenCommit, v *big.Int) (*PedersenCommit, error) {
	commit := new(PedersenCommit)
	commit.Comm = p.Comm.Add(EC.G.Mult(v))
	return commit, nil
}

func PedersenAddCommit(p *PedersenCommit, com *PedersenCommit) (*PedersenCommit, error) {
	commit := new(PedersenCommit)
	commit.Comm = p.Comm.Add(com.Comm)
	return commit, nil
}

func PedersenSubCommit(p *PedersenCommit, com *PedersenCommit) (*PedersenCommit, error) {
	commit := new(PedersenCommit)
	commit.Comm = p.Comm.Add(com.Comm.Neg())
	return commit, nil
}
