package rangproof

import "math/big"

type PedersenCommit struct {
	//G ECPoint // G value for commitments of a single value
	//H ECPoint // H value for commitments of a single value
	Comm ECPoint
}

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
